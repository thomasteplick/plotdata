/*
	Plot data in complex(real, imag) from a file, one space-separated coordinate per line.
	r1 i1
	r2 i2
	...
	rn in

	The html/template package is used to generate the html sent to the client.
	Use CSS display grid to display a 300x300 grid of cells.
	Use CSS flexbox to display the labels on the x and y axes.

	Calculate the power spectral density (PSD) using Welch periodogram spectral
	estimation and plot it.  Average the periodograms using 50% overlap with a
	user-chosen window function such as Bartlett, Hann, Hamming, or Welch.
*/

package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"math/rand"
	"net/http"
	"os"
	"path"
	"strconv"
	"strings"
	"text/template"
	"time"
)

const (
	rows               = 300                                // #rows in grid
	columns            = 300                                // #columns in grid
	tmpltime           = "templates/plottimedata.html"      // html template address
	tmplfrequency      = "templates/plotfrequencydata.html" // html template address
	addr               = "127.0.0.1:8080"                   // http server listen address
	patternPlotData    = "/plotdata"                        // http handler pattern for plotting data
	patternGenerator   = "/generatedata"                    // http handler pattern for data generation
	patternPlotOptions = "/plotoptions"                     // http handler for Plot Options
	xlabels            = 11                                 // # labels on x axis
	ylabels            = 11                                 // # labels on y axis
	dataDir            = "data/"                            // directory for the data files
	deg2rad            = math.Pi / 180.0                    // convert degrees to radians
)

// Type to contain all the HTML template actions
type PlotT struct {
	Grid     []string // plotting grid
	Status   string   // status of the plot
	Xlabel   []string // x-axis labels
	Ylabel   []string // y-axis labels
	Filename string   // filename to plot
}

// Type to hold the minimum and maximum data values
type Endpoints struct {
	xmin float64
	xmax float64
	ymin float64
	ymax float64
}

var (
	timeTmpl *template.Template
	freqTmpl *template.Template
)

// init parses the html template files are done only once
func init() {
	timeTmpl = template.Must(template.ParseFiles(tmpltime))
	freqTmpl = template.Must(template.ParseFiles(tmplfrequency))
}

// findEndpoints finds the minimum and maximum data values
func (ep *Endpoints) findEndpoints(input *bufio.Scanner, rad float64) {
	ep.xmax = -math.MaxFloat64
	ep.xmin = math.MaxFloat64
	ep.ymax = -math.MaxFloat64
	ep.ymin = math.MaxFloat64
	for input.Scan() {
		line := input.Text()
		// Each line has 2 space-separated values: x y coordinates can be int or float
		coord := strings.Split(line, " ")
		var (
			x, y float64
			err  error
		)
		if x, err = strconv.ParseFloat(coord[0], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", coord[0], err)
			continue
		}

		if y, err = strconv.ParseFloat(coord[1], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", coord[1], err)
			continue
		}

		// rotation
		if rad != 0.0 {
			xrot := x*math.Cos(rad) + y*math.Sin(rad)
			yrot := -x*math.Sin(rad) + y*math.Cos(rad)
			x = xrot
			y = yrot
		}

		if x > ep.xmax {
			ep.xmax = x
		}
		if x < ep.xmin {
			ep.xmin = x
		}

		if y > ep.ymax {
			ep.ymax = y
		}
		if y < ep.ymin {
			ep.ymin = y
		}
	}
}

// processTimeDomain plots the time domain data from disk file
func processTimeDomain(w http.ResponseWriter, r *http.Request, filename string) error {

	// main data structure
	var (
		plot      PlotT
		xscale    float64
		yscale    float64
		endpoints Endpoints
	)

	plot.Grid = make([]string, rows*columns)
	plot.Xlabel = make([]string, xlabels)
	plot.Ylabel = make([]string, ylabels)

	// Open file
	f, err := os.Open(filename)
	if err == nil {
		// Mark the data x-y coordinate online at the corresponding
		// grid row/column.
		input := bufio.NewScanner(f)

		// Determine if rotate requested and perform the rotation of x and y
		rotate := r.FormValue("rotate")
		rad := 0.0
		if len(rotate) > 0 {
			deg, err := strconv.ParseFloat(rotate, 64)
			if err != nil {
				plot.Status = "Rotate degree conversion error"
				fmt.Printf("Rotate degree %v conversion error: %v", rotate, err)
			} else {
				rad = deg2rad * deg
			}
		}
		endpoints.findEndpoints(input, rad)

		f.Close()
		f, err = os.Open(filename)
		if err == nil {
			defer f.Close()
			input := bufio.NewScanner(f)

			// Determine if zoom requested and validate endpoints
			zoomxstart := r.FormValue("zoomxstart")
			zoomxend := r.FormValue("zoomxend")
			zoomystart := r.FormValue("zoomystart")
			zoomyend := r.FormValue("zoomyend")
			if len(zoomxstart) > 0 && len(zoomxend) > 0 &&
				len(zoomystart) > 0 && len(zoomyend) > 0 {
				x1, err1 := strconv.ParseFloat(zoomxstart, 64)
				x2, err2 := strconv.ParseFloat(zoomxend, 64)
				y1, err3 := strconv.ParseFloat(zoomystart, 64)
				y2, err4 := strconv.ParseFloat(zoomyend, 64)

				if err1 != nil || err2 != nil || err3 != nil || err4 != nil {
					plot.Status = "Zoom x or y values are not numbers."
					fmt.Printf("Zoom error: x start error = %v, x end error = %v\n", err1, err2)
					fmt.Printf("Zoom error: y start error = %v, y end error = %v\n", err3, err4)
				} else {
					if (x1 < endpoints.xmin || x1 > endpoints.xmax) ||
						(x2 < endpoints.xmin || x2 > endpoints.xmax) || (x1 >= x2) {
						plot.Status = "Zoom values are not in x range."
						fmt.Printf("Zoom error: start or end value not in x range.\n")
					} else if (y1 < endpoints.ymin || y1 > endpoints.ymax) ||
						(y2 < endpoints.ymin || y2 > endpoints.ymax) || (y1 >= y2) {
						plot.Status = "Zoom values are not in y range."
						fmt.Printf("Zoom error: start or end value not in y range.\n")
					} else {
						// Valid Zoom endpoints, replace the previous min and max values
						endpoints.xmin = x1
						endpoints.xmax = x2
						endpoints.ymin = y1
						endpoints.ymax = y2
					}
				}
			}

			// Calculate scale factors for x and y
			xscale = (columns - 1) / (endpoints.xmax - endpoints.xmin)
			yscale = (rows - 1) / (endpoints.ymax - endpoints.ymin)

			for input.Scan() {
				line := input.Text()
				// Each line has 2 space-separated values: x y coordinates can be int or float
				coord := strings.Split(line, " ")
				var x, y float64
				if x, err = strconv.ParseFloat(coord[0], 64); err != nil {
					fmt.Printf("String %s conversion to float error: %v\n", coord[0], err)
					continue
				}

				if y, err = strconv.ParseFloat(coord[1], 64); err != nil {
					fmt.Printf("String %s conversion to float error: %v\n", coord[1], err)
					continue
				}

				// rotation
				if rad != 0.0 {
					xrot := x*math.Cos(rad) + y*math.Sin(rad)
					yrot := -x*math.Sin(rad) + y*math.Cos(rad)
					x = xrot
					y = yrot
				}

				// Check if inside the zoom values
				if x < endpoints.xmin || x > endpoints.xmax || y < endpoints.ymin || y > endpoints.ymax {
					continue
				}

				// This cell location (row,col) is on the line
				//row := (rows - 1) - int((y-endpoints.ymin)*yscale+.5)
				row := int((endpoints.ymax-y)*yscale + .5)
				col := int((x-endpoints.xmin)*xscale + .5)
				plot.Grid[row*columns+col] = "online"
			}
			// Set plot status if no errors
			if len(plot.Status) == 0 {
				plot.Status = fmt.Sprintf("Status: Data plotted from (%.3f,%.3f) to (%.3f,%.3f)",
					endpoints.xmin, endpoints.ymin, endpoints.xmax, endpoints.ymax)
			}

		} else {
			// Set plot status
			fmt.Printf("Error opening file %s: %v\n", filename, err)
			return fmt.Errorf("error opening file %s: %v", filename, err)
		}
	} else {
		// Set plot status
		fmt.Printf("Error opening file %s: %v\n", filename, err)
		return fmt.Errorf("error opening file %s: %v", filename, err)
	}

	// Construct x-axis labels
	incr := (endpoints.xmax - endpoints.xmin) / (xlabels - 1)
	x := endpoints.xmin
	// First label is empty for alignment purposes
	for i := range plot.Xlabel {
		plot.Xlabel[i] = fmt.Sprintf("%.2f", x)
		x += incr
	}

	// Construct the y-axis labels
	incr = (endpoints.ymax - endpoints.ymin) / (ylabels - 1)
	y := endpoints.ymin
	for i := range plot.Ylabel {
		plot.Ylabel[i] = fmt.Sprintf("%.2f", y)
		y += incr
	}

	// Enter the filename in the form
	plot.Filename = path.Base(filename)

	// Write to HTTP using template and grid
	if err := timeTmpl.Execute(w, plot); err != nil {
		log.Fatalf("Write to HTTP output using template with grid error: %v\n", err)
	}

	return nil
}

// processFrequencyDomain calculates the power spectral density (PSD) and plots it
func processFrequencyDomain(w http.ResponseWriter, r *http.Request, filename string) error {
	// Use complex128 for FFT computation
	// Get the number of complex samples nn, open file and count lines, close the file
	// Get FFT size and number of segments K in PSD periodogram from form
	// Require 1 <= K <= 10
	// Get window type:  Bartlett, Welch, Hann, Hamming, etc
	// Segment size m = nn/(K+1)
	// Check FFT Size >= 2*m, 50%  overlap of segments
	// Check FFT Size is a power of 2:  2^n => n = int(math.Log2(float64(FFT Size) + .5)), this rounds up to
	//     nearest FFT size that is a power of 2.
	// FFT Size = math.Exp2(float64(n)), FFT Size = math.Max(2*m, FFT Size)

	// Reopen the data file
	// Read in m samples into buf[m]
	// loop over file samples
	//     copy m samples from buf[m] to front of buf[2m]
	//     read in m file samples into buf[m]
	//     copy m samples from buf[m] to back half of buf[2m]
	//     window the 2m samples with chosen window
	//     zero-pad N-2m in buf[N], float64
	//     Perform N-point complex FFT and add squares to previous values in PSD[N/2] -> buf[n], buf[N-n]
	// Normalize the PSD using N*Sum(w[i]*w[i])*K
	// Store in the form:  FFT Size, window type, number of samples nn, K segments
	// Plot the PSD N/2 float64 values, execute the data on the plotfrequency.html template

	return nil
}

// handlePlotData opens the file, does frequency domain processing if required,
// and plots the data in the time or frequency domain
func handlePlotData(w http.ResponseWriter, r *http.Request) {
	// main data structure
	var (
		plot      PlotT
		endpoints Endpoints
	)
	plot.Grid = make([]string, rows*columns)
	plot.Xlabel = make([]string, xlabels)
	plot.Ylabel = make([]string, ylabels)

	filename := r.FormValue("filename")
	// choose time or frequency domain processing
	if len(filename) > 0 {

		domain := r.FormValue("domain")
		switch domain {
		case "time":
			err := processTimeDomain(w, r, "data/"+filename)
			if err != nil {
				plot.Status = err.Error()
			}
		case "frequency":
			err := processFrequencyDomain(w, r, "data/"+filename)
			if err != nil {
				plot.Status = err.Error()
			}
		default:
			plot.Status = fmt.Sprintf("Invalid domain choice: %s", domain)
		}
	} else {
		plot.Status = "Status: Provide a data file to plot."
		// Construct x-axis labels
		incr := (endpoints.xmax - endpoints.xmin) / (xlabels - 1)
		x := endpoints.xmin
		// First label is empty for alignment purposes
		for i := range plot.Xlabel {
			plot.Xlabel[i] = fmt.Sprintf("%.2f", x)
			x += incr
		}

		// Construct the y-axis labels
		incr = (endpoints.ymax - endpoints.ymin) / (ylabels - 1)
		y := endpoints.ymin
		for i := range plot.Ylabel {
			plot.Ylabel[i] = fmt.Sprintf("%.2f", y)
			y += incr
		}

		// Write to HTTP using template and grid
		if err := timeTmpl.Execute(w, plot); err != nil {
			log.Fatalf("Write to HTTP output using template with grid error: %v\n", err)
		}
	}
}

// HTTP handler for /plotoptions connections
func handlePlotOptions(w http.ResponseWriter, r *http.Request) {
	http.ServeFile(w, r, "templates/plotoptions.html")
}

// handleGenerating creates x-y data and saves to disk files
func handleGenerating(w http.ResponseWriter, r *http.Request) {
	const (
		samples         = 800
		cycles          = 4
		k       float64 = cycles * math.Pi
	)

	// 1. sine wave
	f, err := os.Create(path.Join(dataDir, "sin.txt"))
	if err != nil {
		log.Fatalf("Create %v error: %v", path.Join(dataDir, "sin.txt"), err)
	}

	xstart := -1.0
	xend := 1.0
	xincr := (xend - xstart) / float64(samples)
	for x := xstart; x < xend; x += xincr {
		fmt.Fprintf(f, "%v %v\n", x, math.Sin(k*x))
	}
	fmt.Fprintf(w, "Wrote %v samples to %v\n", samples, path.Join(dataDir, "sin.txt"))
	f.Close()

	// 2. circle
	f, err = os.Create(path.Join(dataDir, "circle.txt"))
	if err != nil {
		log.Fatalf("Create %v error: %v", path.Join(dataDir, "circle.txt"), err)
	}
	// radius
	rad := 8.0
	incr := 2.0 * rad / samples
	x := -rad
	var y float64
	for i := 0; i < samples; i++ {
		y = math.Sqrt(rad*rad - x*x)
		fmt.Fprintf(f, "%v %v\n", x, y)
		fmt.Fprintf(f, "%v %v\n", x, -y)
		x += incr
	}
	fmt.Fprintf(w, "Wrote %v samples to %v\n", samples, path.Join(dataDir, "circle.txt"))
	f.Close()

	// 3. decaying sawtooth
	f, err = os.Create(path.Join(dataDir, "sawtooth.txt"))
	if err != nil {
		log.Fatalf("Create %v error: %v", path.Join(dataDir, "sawtooth.txt"), err)
	}
	y = 100.0
	x = 25.0
	xincr = 1.5
	yincr := 3.0 * xincr
	for c := 0; c < cycles; c++ {
		samplesPerCycle := samples / cycles
		for n := 0; n < samplesPerCycle; n++ {
			fmt.Fprintf(f, "%v %v\n", x, y)
			// reverse y direction
			if n == samplesPerCycle/2 {
				yincr = -yincr
			}
			x += xincr
			y += yincr
		}
		// decay and reverse direction
		yincr *= -.8
	}
	fmt.Fprintf(w, "Wrote %v samples to %v\n", samples, path.Join(dataDir, "sawtooth.txt"))
	f.Close()

	// 4. Random data
	f, err = os.Create(path.Join(dataDir, "random.txt"))
	if err != nil {
		log.Fatalf("Create %v error: %v", path.Join(dataDir, "random.txt"), err)
	}
	rand.Seed(time.Now().Unix())
	m := 1000.0
	for i := 0; i < samples; i++ {
		x := m * rand.Float64()
		y := m * rand.Float64()
		fmt.Fprintf(f, "%v %v\n", x, y)
	}
	fmt.Fprintf(w, "Wrote %v samples to %v\n", samples, path.Join(dataDir, "random.txt"))
	f.Close()

	// 5. Spiral
	f, err = os.Create(path.Join(dataDir, "spiral.txt"))
	if err != nil {
		log.Fatalf("Create %v error: %v", path.Join(dataDir, "spiral.txt"), err)
	}

	rad = 0.0
	ang := 0.0
	radIncr := .1
	angIncr := 2.0 * math.Pi / (samples / 2.0)
	for s := 0; s < cycles*samples; s++ {
		x := rad * math.Cos(ang)
		y := rad * math.Sin(ang)
		fmt.Fprintf(f, "%v %v\n", x, y)
		rad += radIncr
		ang += angIncr
	}
	fmt.Fprintf(w, "Wrote %v samples to %v\n", cycles*samples, path.Join(dataDir, "spiral.txt"))
	f.Close()

	// 6. Cardoid r = K(1-cos(ang)), where K is a positive constant
	f, err = os.Create(path.Join(dataDir, "cardoid.txt"))
	if err != nil {
		log.Fatalf("Create %v error: %v", path.Join(dataDir, "cardoid.txt"), err)
	}

	const K = 13.0
	angIncr = 2.0 * math.Pi / samples
	for ang := 0.0; ang < 2*math.Pi; ang += angIncr {
		r := K * (1.0 - math.Cos(ang))
		x := r * math.Cos(ang)
		y := r * math.Sin(ang)
		fmt.Fprintf(f, "%v %v\n", x, y)
	}
	fmt.Fprintf(w, "Wrote %v samples to %v\n", samples, path.Join(dataDir, "cardoid.txt"))
	f.Close()
}

// executive program
func main() {
	// Setup http server with handler for reading file time data, calculating frequency domain data,
	// and plotting time or frequency domain data
	http.HandleFunc(patternPlotData, handlePlotData)

	// Setup http server with handler for generating data for testing
	http.HandleFunc(patternGenerator, handleGenerating)

	// Setup http server with handler for choosing data file and either time or frequency domain plot
	http.HandleFunc(patternPlotOptions, handlePlotOptions)
	fmt.Printf("Plot Data Server listening on %v.\n", addr)

	http.ListenAndServe(addr, nil)
}
