/*
	Plot data in complex(real, imag) from a file, one complex number  per line.
	r1 i1
	r2 i2
	...
	rn in

	The html/template package is used to generate the html sent to the client.
	Use CSS display grid to display a 300x300 grid of cells.
	Use CSS flexbox to display the labels on the x and y axes.

	Calculate the power spectral density (PSD) using Welch periodogram spectral
	estimation and plot it.  Average the periodograms using 50% overlap with a
	user-chosen window function such as Bartlett, Hanning, Hamming, or Welch.
*/

package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"math/rand"
	"net/http"
	"os"
	"path"
	"strconv"
	"strings"
	"text/template"

	"github.com/mjibson/go-dsp/fft"
)

const (
	rows                = 300                                // #rows in grid
	columns             = 300                                // #columns in grid
	tmpltime            = "templates/plottimedata.html"      // html template address
	tmplfrequency       = "templates/plotfrequencydata.html" // html template address
	tmploptions         = "templates/plotoptions.html"       // html template address
	tmplgenerate        = "templates/generatedata.html"      // html template address
	addr                = "127.0.0.1:8080"                   // http server listen address
	patternPlotData     = "/plotdata"                        // http handler pattern for plotting data
	patternGenerateData = "/generatedata"                    // http handler pattern for data generation
	patternPlotOptions  = "/plotoptions"                     // http handler for Plot Options
	xlabels             = 11                                 // # labels on x axis
	ylabels             = 11                                 // # labels on y axis
	dataDir             = "data/"                            // directory for the data files
	deg2rad             = math.Pi / 180.0                    // convert degrees to radians
)

// Type to contain all the HTML template actions
type PlotT struct {
	Grid         []string // plotting grid
	Status       string   // status of the plot
	Xlabel       []string // x-axis labels
	Ylabel       []string // y-axis labels
	Filename     string   // filename to plot
	SamplingFreq string   // data sampling rate in Hz
	FFTSegments  string   // FFT segments, K
	FFTSize      string   // FFT size
	Samples      string   // complex samples in data file
}

// Type to hold the minimum and maximum data values
type Endpoints struct {
	xmin float64
	xmax float64
	ymin float64
	ymax float64
}

// Window function type
type Window func(n int, m int) complex128

var (
	timeTmpl     *template.Template
	freqTmpl     *template.Template
	generateTmpl *template.Template
	optionTmpl   *template.Template
	winType      = []string{"Bartlett", "Welch", "Hamming", "Hanning"}
)

// Bartlett window
func bartlett(n int, m int) complex128 {
	real := 1.0 - math.Abs((float64(n)-float64(m))/float64(m))
	return complex(real, 0)
}

// Welch window
func welch(n int, m int) complex128 {
	x := math.Abs((float64(n) - float64(m)) / float64(m))
	real := 1.0 - x*x
	return complex(real, 0)
}

// Hamming window
func hamming(n int, m int) complex128 {
	return complex(.54-.46*math.Cos(math.Pi*float64(n)/float64(m)), 0)
}

// Hanning window
func hanning(n int, m int) complex128 {
	return complex(.5-.5*math.Cos(math.Pi*float64(n)/float64(m)), 0)
}

// Rectangle window
func rectangle(n int, m int) complex128 {
	return 1.0
}

// init parses the html template files are done only once
func init() {
	timeTmpl = template.Must(template.ParseFiles(tmpltime))
	freqTmpl = template.Must(template.ParseFiles(tmplfrequency))
	optionTmpl = template.Must(template.ParseFiles(tmploptions))
	generateTmpl = template.Must(template.ParseFiles(tmplgenerate))
}

// findEndpoints finds the minimum and maximum data values
func (ep *Endpoints) findEndpoints(input *bufio.Scanner, rad float64) {
	ep.xmax = -math.MaxFloat64
	ep.xmin = math.MaxFloat64
	ep.ymax = -math.MaxFloat64
	ep.ymin = math.MaxFloat64
	for input.Scan() {
		line := input.Text()
		// Each line has 2 or 3 space-separated values, depending on if it is real or complex data:
		// time real, or x y for euclidean data
		// time real imaginary for complex data
		values := strings.Split(line, " ")
		var (
			x, y, t float64
			err     error
		)
		// real data
		if len(values) == 2 {
			if x, err = strconv.ParseFloat(values[0], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
				continue
			}

			if y, err = strconv.ParseFloat(values[1], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
				continue
			}
			// complex data
		} else if len(values) == 3 {
			if t, err = strconv.ParseFloat(values[0], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
				continue
			}

			if x, err = strconv.ParseFloat(values[1], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
				continue
			}

			if y, err = strconv.ParseFloat(values[2], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[2], err)
				continue
			}

			// Calculate the modulus of the complex data becomes the y-axis
			// The time becomes the x-axis
			y = math.Sqrt(x*x + y*y)
			x = t
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

// gridFill inserts the data points in the grid
func gridFill(plot *PlotT, xscale float64, yscale float64, endpoints Endpoints, rad float64, input *bufio.Scanner) error {
	for input.Scan() {
		line := input.Text()
		// Each line has 2 or 3 space-separated values, depending on if it is real or complex data:
		// time real, or x y for euclidean data
		// time real imaginary for complex data
		values := strings.Split(line, " ")
		var (
			x, y, t float64
			err     error
		)
		// real data
		if len(values) == 2 {
			if x, err = strconv.ParseFloat(values[0], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
				return err
			}

			if y, err = strconv.ParseFloat(values[1], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
				return err
			}
			// complex data
		} else if len(values) == 3 {
			if t, err = strconv.ParseFloat(values[0], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
				return err
			}

			if x, err = strconv.ParseFloat(values[1], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
				return err
			}

			if y, err = strconv.ParseFloat(values[2], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[2], err)
				return err
			}

			// Calculate the modulus of the complex data becomes the y-axis
			// The time becomes the x-axis
			y = math.Sqrt(x*x + y*y)
			x = t
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
		row := int((endpoints.ymax-y)*yscale + .5)
		col := int((x-endpoints.xmin)*xscale + .5)
		plot.Grid[row*columns+col] = "online"
	}
	return nil
}

// gridFillInterp inserts the data points in the grid and draws a straight line between points
func gridFillInterp(plot *PlotT, xscale float64, yscale float64, endpoints Endpoints, rad float64, input *bufio.Scanner) error {

	var (
		x, y, t      float64
		prevX, prevY float64
		prevSet      bool = true
		err          error
	)

	const lessen = 10

	// Get first sample
	input.Scan()
	line := input.Text()
	// Each line has 2 or 3 space-separated values, depending on if it is real or complex data:
	// time real, or x y for euclidean data
	// time real imaginary for complex data
	values := strings.Split(line, " ")
	// real data
	if len(values) == 2 {
		if x, err = strconv.ParseFloat(values[0], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
			return err
		}

		if y, err = strconv.ParseFloat(values[1], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
			return err
		}
		// complex data
	} else if len(values) == 3 {
		if t, err = strconv.ParseFloat(values[0], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
			return err
		}

		if x, err = strconv.ParseFloat(values[1], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
			return err
		}

		if y, err = strconv.ParseFloat(values[2], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[2], err)
			return err
		}

		// Calculate the modulus of the complex data becomes the y-axis
		// The time becomes the x-axis
		y = math.Sqrt(x*x + y*y)
		x = t
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
		prevSet = false
	} else {
		// This cell location (row,col) is on the line
		row := int((endpoints.ymax-y)*yscale + .5)
		col := int((x-endpoints.xmin)*xscale + .5)
		plot.Grid[row*columns+col] = "online"

		prevX = x
		prevY = y
	}

	// Scale factor to determine the number of interpolation points
	lenEP := math.Sqrt((endpoints.xmax-endpoints.xmin)*(endpoints.xmax-endpoints.xmin) +
		(endpoints.ymax-endpoints.ymin)*(endpoints.ymax-endpoints.ymin))

	// Continue with the rest of the points in the file
	for input.Scan() {
		line = input.Text()
		// Each line has 2 or 3 space-separated values, depending on if it is real or complex data:
		// time real, or x y for euclidean data
		// time real imaginary for complex data
		values := strings.Split(line, " ")
		// real data
		if len(values) == 2 {
			if x, err = strconv.ParseFloat(values[0], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
				return err
			}

			if y, err = strconv.ParseFloat(values[1], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
				return err
			}
			// complex data
		} else if len(values) == 3 {
			if t, err = strconv.ParseFloat(values[0], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
				return err
			}

			if x, err = strconv.ParseFloat(values[1], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
				return err
			}

			if y, err = strconv.ParseFloat(values[2], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[2], err)
				return err
			}

			// Calculate the modulus of the complex data becomes the y-axis
			// The time becomes the x-axis
			y = math.Sqrt(x*x + y*y)
			x = t
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
		} else if !prevSet {
			prevSet = true
			prevX = x
			prevY = y
		}

		// This cell location (row,col) is on the line
		row := int((endpoints.ymax-y)*yscale + .5)
		col := int((x-endpoints.xmin)*xscale + .5)
		plot.Grid[row*columns+col] = "online"

		// Interpolate the points between previous point and current point

		lenEdge := math.Sqrt((x-prevX)*(x-prevX) + (y-prevY)*(y-prevY))
		ncells := int(columns*lenEdge/lenEP) / lessen // number of points to interpolate
		stepX := (x - prevX) / float64(ncells)
		stepY := (y - prevY) / float64(ncells)

		// loop to draw the points
		interpX := prevX
		interpY := prevY
		for i := 0; i < ncells; i++ {
			row := int((endpoints.ymax-interpY)*yscale + .5)
			col := int((interpX-endpoints.xmin)*xscale + .5)
			plot.Grid[row*columns+col] = "online"
			interpX += stepX
			interpY += stepY
		}

		// Update the previous point with the current point
		prevX = x
		prevY = y
	}
	return nil
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

			// Check for interpolation and fill in the grid with the data points
			interp := r.FormValue("interpolate")
			if interp == "interpolate" {
				err = gridFillInterp(&plot, xscale, yscale, endpoints, rad, input)
			} else {
				err = gridFill(&plot, xscale, yscale, endpoints, rad, input)
			}
			if err != nil {
				return err
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
		log.Fatalf("Write to HTTP output using template with error: %v\n", err)
	}

	return nil
}

// processFrequencyDomain calculates the power spectral density (PSD) and plots it
func processFrequencyDomain(w http.ResponseWriter, r *http.Request, filename string) error {
	// Use complex128 for FFT computation
	// Get the number of complex samples nn, open file and count lines, close the file

	var (
		plot          PlotT // main data structure to execute with parsed html template
		endpoints     Endpoints
		N             int                                                        //  complex FFT size
		nn            int                                                        // number of complex samples in the data file
		K             int                                                        //  number of segments used in PSD with 50% overlap
		m             int                                                        // complex segment size
		win           string                                                     // FFT window type
		window        = make(map[string]Window, len(winType))                    // map of window functions
		sumWindow     float64                                                    // sum of squared window values for normalization
		normalizerPSD float64                                                    // normalizer for PSD
		PSD           []float64                                                  // power spectral density
		psdMax        float64                                 = -math.MaxFloat64 // maximum PSD value
		psdMin        float64                                 = math.MaxFloat64  // minimum PSD value
		xscale        float64                                                    // data to grid in x direction
		yscale        float64                                                    // data to grid in y direction
		samplingRate  float64                                                    // sampling rate in Hz
	)

	plot.Grid = make([]string, rows*columns)
	plot.Xlabel = make([]string, xlabels)
	plot.Ylabel = make([]string, ylabels)

	// Put the window functions in the map
	window["Bartlett"] = bartlett
	window["Welch"] = welch
	window["Hamming"] = hamming
	window["Hanning"] = hanning
	window["Rectangle"] = rectangle

	// Open file
	f, err := os.Open(filename)
	if err == nil {
		input := bufio.NewScanner(f)
		// Number of real or complex data samples
		for input.Scan() {
			line := input.Text()
			if len(line) > 0 {
				nn++
			}
		}
		fmt.Printf("Data file %s has %d samples\n", filename, nn)

		f.Close()

		// Get number of segments from HTML form
		// Number of segments to average the periodograms to reduce the variance
		tmp := r.FormValue("fftsegments")
		if len(tmp) == 0 {
			return fmt.Errorf("enter number of FFT segments")
		}
		K, err = strconv.Atoi(tmp)
		if err != nil {
			fmt.Printf("FFT segments string convert error: %v\n", err)
			return fmt.Errorf("fft segments string convert error: %v", err)
		}

		// Require 1 <= K <= 20
		if K < 1 {
			K = 1
		} else if K > 20 {
			K = 20
		}

		// segment size complex samples
		m = nn / (K + 1)

		// Get window type:  Bartlett, Welch, Hanning, Hamming, etc
		// Multiply the samples by the window to reduce spectral leakage
		// caused by high sidelobes in rectangular window
		win = r.FormValue("fftwindow")
		if len(win) == 0 {
			return fmt.Errorf("enter FFT window type")
		}
		w, ok := window[win]
		if !ok {
			fmt.Printf("Invalid FFT window type: %v\n", win)
			return fmt.Errorf("invalid FFT window type: %v", win)
		}
		// sum the window values for PSD normalization due to windowing
		for i := 0; i < 2*m; i++ {
			x := cmplx.Abs(w(i, m))
			sumWindow += x * x
		}
		fmt.Printf("%s window sum = %.2f\n", win, sumWindow)

		// Get FFT size from HTML form
		// Check FFT Size >= 2*m, using 50%  overlap of segments
		// Check FFT Size is a power of 2:  2^n
		tmp = r.FormValue("fftsize")
		if len(tmp) == 0 {
			return fmt.Errorf("enter FFT size")
		}
		N, err = strconv.Atoi(tmp)
		if err != nil {
			return fmt.Errorf("fft size string convert error: %v", err)
		}

		if N < rows {
			fmt.Printf("FFT size < %d\n", rows)
			N = rows
		} else if N > rows*rows {
			fmt.Printf("FFT size > %d\n", rows*rows)
			N = rows * rows
		}
		// This rounds up to nearest FFT size that is a power of 2
		N = int(math.Exp2(float64(int(math.Log2(float64(N)) + .5))))
		fmt.Printf("N=%v\n", N)

		if N < 2*m {
			fmt.Printf("FFT Size %d not greater than 2%d\n", N, 2*m)
			return fmt.Errorf("fft Size %d not greater than 2*%d", N, 2*m)
		}

		// Power Spectral Density, PSD[N/2] is the Nyquist critical frequency
		// It is the (sampling frequency)/2, the highest non-aliased frequency
		PSD = make([]float64, N/2+1)

		// Reopen the data file
		f, err = os.Open(filename)
		if err == nil {
			defer f.Close()
			bufm := make([]complex128, m)
			bufN := make([]complex128, N)
			input := bufio.NewScanner(f)
			// Read in initial m samples into buf[m] to start the processing loop
			diff := 0.0
			prev := 0.0
			for k := 0; k < m; k++ {
				input.Scan()
				line := input.Text()
				// Each line has 2 or 3 space-separated values, depending on if it is real or complex data:
				// time real
				// time real imaginary for complex data
				values := strings.Split(line, " ")
				// real data
				if len(values) == 2 {
					// time real, calculate the sampling rate from the time steps
					var t, r float64

					if t, err = strconv.ParseFloat(values[0], 64); err != nil {
						fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
						continue
					}
					if r, err = strconv.ParseFloat(values[1], 64); err != nil {
						fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
						continue
					}

					if k == 0 {
						prev = t
					} else {
						diff += t - prev
						prev = t
					}
					bufm[k] = complex(r, 0)

					// complex data
				} else if len(values) == 3 {
					// time real imag
					var t, r, i float64

					// calculate the sampling rate from the time steps
					if t, err = strconv.ParseFloat(values[0], 64); err != nil {
						fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
						continue
					}
					if r, err = strconv.ParseFloat(values[1], 64); err != nil {
						fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
						continue
					}

					if i, err = strconv.ParseFloat(values[2], 64); err != nil {
						fmt.Printf("String %s conversion to float error: %v\n", values[2], err)
						continue
					}

					if k == 0 {
						prev = t
					} else {
						diff += t - prev
						prev = t
					}
					bufm[k] = complex(r, i)
				}
			}

			// Average the time steps and invert to get the sampling rate
			samplingRate = 1.0 / (diff / float64(m-1))
			fmt.Printf("sampling rate = %.2f\n", samplingRate)

			scanOK := true
			// loop over the rest of the file, reading in m samples at a time until EOF
			for {
				// Put the previous m samples in the front of the buffer N to
				// overlap segments
				copy(bufN, bufm)
				// Put the next m samples in back of the previous m samples
				kk := 0
				for k := 0; k < m; k++ {
					scanOK = input.Scan()
					if !scanOK {
						break
					}
					line := input.Text()
					// Each line has 2 space-separated values: real and imaginary
					values := strings.Split(line, " ")
					// real data
					if len(values) == 2 {
						// time real, don't need the time
						var r float64

						if r, err = strconv.ParseFloat(values[1], 64); err != nil {
							fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
							continue
						}
						// real data, imaginary component is zero
						bufm[k] = complex(r, 0)

						// complex data
					} else if len(values) == 3 {
						// time real imag
						var r, i float64
						// Don't need the time
						if r, err = strconv.ParseFloat(values[1], 64); err != nil {
							fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
							continue
						}

						if i, err = strconv.ParseFloat(values[2], 64); err != nil {
							fmt.Printf("String %s conversion to float error: %v\n", values[2], err)
							continue
						}

						bufm[k] = complex(r, i)
					}
					kk++
				}
				// Check for the normal EOF and the abnormal scan error
				// EOF does not give an error but is considered normal termination
				if !scanOK {
					if input.Err() != nil {
						fmt.Printf("Data file scan error: %v", input.Err().Error())
						return fmt.Errorf("data file scan error: %v", input.Err())
					}
				}
				// Put the next kk samples in back of the previous m
				copy(bufN[m:], bufm[:kk])

				// window the m + kk samples with chosen window
				for i := 0; i < m+kk; i++ {
					bufN[i] *= w(i, m)
				}

				// zero-pad N-m-kk samples in buf[N]
				for i := m + kk; i < N; i++ {
					bufN[i] = 0
				}

				// Perform N-point complex FFT and add squares to previous values in PSD[N/2+1]
				fourierN := fft.FFT(bufN)
				x := cmplx.Abs(fourierN[0])
				PSD[0] += x * x
				for i := 1; i < N/2; i++ {
					// Use positive and negative frequencies -> bufN[N-i] = bufN[-i]
					xi := cmplx.Abs(fourierN[i])
					xNi := cmplx.Abs(fourierN[N-i])
					PSD[i] += xi*xi + xNi*xNi
				}
				x = cmplx.Abs(fourierN[N/2])
				PSD[N/2] += x * x

				// part of K*Sum(w[i]*w[i]) PSD normalizer
				normalizerPSD += sumWindow

				// EOF reached
				if !scanOK {
					break
				}
			} // K segments done

			// Normalize the PSD using K*Sum(w[i]*w[i])
			// Use log plot for wide dynamic range
			if r.FormValue("plottype") == "linear" {
				for i := range PSD {
					PSD[i] /= normalizerPSD
					if PSD[i] > psdMax {
						psdMax = PSD[i]
					}
					if PSD[i] < psdMin {
						psdMin = PSD[i]
					}
				}
				// log10 in dB
			} else {
				for i := range PSD {
					PSD[i] /= normalizerPSD
					PSD[i] = 10.0 * math.Log10(PSD[i])
					if PSD[i] > psdMax {
						psdMax = PSD[i]
					}
					if PSD[i] < psdMin {
						psdMin = PSD[i]
					}
				}
			}

			endpoints.xmin = 0
			endpoints.xmax = float64(N / 2) // equivalent to Nyquist critical frequency
			endpoints.ymin = psdMin
			endpoints.ymax = psdMax

			// Calculate scale factors for x and y
			xscale = (columns - 1) / (endpoints.xmax - endpoints.xmin)
			yscale = (rows - 1) / (endpoints.ymax - endpoints.ymin)

			// Store the PSD in the plot Grid
			for bin, pow := range PSD {
				row := int((endpoints.ymax-pow)*yscale + .5)
				col := int((float64(bin)-endpoints.xmin)*xscale + .5)
				plot.Grid[row*columns+col] = "online"
			}

			// Store in the form:  FFT Size, window type, number of samples nn, K segments, sampling frequency
			// Plot the PSD N/2 float64 values, execute the data on the plotfrequency.html template

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

	// Apply the  sampling rate in Hz to the x-axis using a scale factor
	// Convert the fft size to sampling rate/2, the Nyquist critical frequency
	sf := 0.5 * samplingRate / endpoints.xmax

	// Construct x-axis labels
	incr := (endpoints.xmax - endpoints.xmin) / (xlabels - 1)
	x := endpoints.xmin
	// First label is empty for alignment purposes
	for i := range plot.Xlabel {
		plot.Xlabel[i] = fmt.Sprintf("%.0f", x*sf)
		x += incr
	}

	// Construct the y-axis labels
	incr = (endpoints.ymax - endpoints.ymin) / (ylabels - 1)
	y := endpoints.ymin
	for i := range plot.Ylabel {
		plot.Ylabel[i] = fmt.Sprintf("%.2f", y)
		y += incr
	}

	// Insert frequency domain parameters in the form
	plot.SamplingFreq = fmt.Sprintf("%.0f", samplingRate)
	plot.FFTSegments = strconv.Itoa(K)
	plot.FFTSize = strconv.Itoa(N)
	plot.Samples = strconv.Itoa(nn)

	// Enter the filename in the form
	plot.Filename = path.Base(filename)

	// Write to HTTP using template and grid
	if err := freqTmpl.Execute(w, plot); err != nil {
		log.Fatalf("Write to HTTP output using template with error: %v\n", err)
	}

	return nil
}

// handlePlotData opens the file, does frequency domain processing if required,
// and plots the data in the time or frequency domain
func handlePlotData(w http.ResponseWriter, r *http.Request) {
	// main data structure
	var (
		plot PlotT
		err  error = nil
	)

	filename := r.FormValue("filename")
	// choose time or frequency domain processing
	if len(filename) > 0 {

		domain := r.FormValue("domain")
		switch domain {
		case "time":
			err = processTimeDomain(w, r, path.Join(dataDir, filename))
			if err != nil {
				plot.Status = err.Error()
			}
		case "frequency":
			err = processFrequencyDomain(w, r, path.Join(dataDir, filename))
			if err != nil {
				plot.Status = err.Error()
			}
		default:
			plot.Status = fmt.Sprintf("Invalid domain choice: %s", domain)
		}

		if err != nil {

			// Write to HTTP using template and grid``
			if err := optionTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
		}
	} else {
		plot.Status = "Missing filename to plot"
		// Write to HTTP using template and grid``
		if err := optionTmpl.Execute(w, plot); err != nil {
			log.Fatalf("Write to HTTP output using template with error: %v\n", err)
		}
	}
}

// HTTP handler for /plotoptions connections
func handlePlotOptions(w http.ResponseWriter, r *http.Request) {
	var plot PlotT
	plot.Status = "Enter filename and plot options"
	// Write to HTTP using template and grid``
	if err := optionTmpl.Execute(w, plot); err != nil {
		log.Fatalf("Write to HTTP output using template with error: %v\n", err)
	}
}

// processFigureData generates and saves to disk the desired figures
func processFigureData(w http.ResponseWriter, r *http.Request, samples int, samplerate int) error {

	// 1. sine wave
	if figure := r.FormValue("figuresine"); figure == "sine" {
		f, err := os.Create(path.Join(dataDir, "sin.txt"))
		if err != nil {
			return fmt.Errorf("create %v error: %v", path.Join(dataDir, "sin.txt"), err)
		}

		// Put the sinusoid at a low frequency so that is discernible when plotted in time domain
		freq := float64(samplerate) / 1000.0
		omega := 2.0 * math.Pi * freq
		delta := 1.0 / float64(samplerate)
		for t := 0.0; t < float64(samples)*delta; t += delta {
			fmt.Fprintf(f, "%v %v\n", t, math.Sin(omega*t))
		}
		f.Close()
	}

	// 2. circle
	if figure := r.FormValue("figurecircle"); figure == "circle" {
		f, err := os.Create(path.Join(dataDir, "circle.txt"))
		if err != nil {
			return fmt.Errorf("create %v error: %v", path.Join(dataDir, "circle.txt"), err)
		}
		// radius
		rad := 8.0
		incr := 2.0 * rad / float64(samples)
		x := -rad
		var y float64
		for i := 0; i < samples; i++ {
			y = math.Sqrt(rad*rad - x*x)
			fmt.Fprintf(f, "%v %v\n", x, y)
			fmt.Fprintf(f, "%v %v\n", x, -y)
			x += incr
		}
		f.Close()
	}

	// 3. decaying sawtooth
	if figure := r.FormValue("figuresawtooth"); figure == "sawtooth" {
		cycles := 4
		f, err := os.Create(path.Join(dataDir, "sawtooth.txt"))
		if err != nil {
			return fmt.Errorf("create %v error: %v", path.Join(dataDir, "sawtooth.txt"), err)
		}
		y := 100.0
		x := 25.0
		xincr := 1.5
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
		f.Close()
	}

	// 4. Spiral
	if figure := r.FormValue("figurespiral"); figure == "spiral" {
		cycles := 4
		f, err := os.Create(path.Join(dataDir, "spiral.txt"))
		if err != nil {
			return fmt.Errorf("create %v error: %v", path.Join(dataDir, "spiral.txt"), err)
		}

		rad := 0.0
		ang := 0.0
		radIncr := .1
		angIncr := 2.0 * math.Pi / (float64(samples) / 2.0)
		for s := 0; s < cycles*samples; s++ {
			x := rad * math.Cos(ang)
			y := rad * math.Sin(ang)
			fmt.Fprintf(f, "%v %v\n", x, y)
			rad += radIncr
			ang += angIncr
		}
		f.Close()
	}

	// 5. Cardoid r = K(1-cos(ang)), where K is a positive constant
	if figure := r.FormValue("figurecardoid"); figure == "cardoid" {
		f, err := os.Create(path.Join(dataDir, "cardoid.txt"))
		if err != nil {
			return fmt.Errorf("create %v error: %v", path.Join(dataDir, "cardoid.txt"), err)
		}

		const K = 13.0
		angIncr := 2.0 * math.Pi / float64(samples)
		for ang := 0.0; ang < 2*math.Pi; ang += angIncr {
			r := K * (1.0 - math.Cos(ang))
			x := r * math.Cos(ang)
			y := r * math.Sin(ang)
			fmt.Fprintf(f, "%v %v\n", x, y)
		}
		f.Close()
	}

	return nil
}

// generateSinsuoids creates a sum of sine waves in Gaussian noise at SNR
func generateSinsuoids(w http.ResponseWriter, r *http.Request, samples int, samplerate int) error {
	// get snr
	temp := r.FormValue("SSsnr")
	if len(temp) == 0 {
		return fmt.Errorf("missing SNR for Sum of Sinsuoids")
	}
	snr, err := strconv.Atoi(temp)
	if err != nil {
		return err
	}

	type Sine struct {
		freq int
		amp  int
	}

	var (
		maxampl  int      = 0 // maximum sine amplitude
		noiseSD  float64      // noise standard deviation
		freqName []string = []string{"SSfreq1", "SSfreq2",
			"SSfreq3", "SSfreq4", "SSfreq5"}
		ampName []string = []string{"SSamp1", "SSamp2", "SSamp3",
			"SSamp4", "SSamp5"}
		signal []Sine = make([]Sine, 0) // sines to generate
	)

	// get the sine frequencies and amplitudes, 1 to 5 possible
	for i, name := range freqName {
		a := r.FormValue(ampName[i])
		f := r.FormValue(name)
		if len(a) > 0 && len(f) > 0 {
			freq, err := strconv.Atoi(f)
			if err != nil {
				return err
			}
			amp, err := strconv.Atoi(a)
			if err != nil {
				return err
			}
			signal = append(signal, Sine{freq: freq, amp: amp})
			if amp > maxampl {
				maxampl = amp
			}
		}
	}
	// Require at least one sine to create
	if len(signal) == 0 {
		return fmt.Errorf("enter frequency and amplitude of 1 to 5 sinsuoids")
	}

	// Calculate the noise standard deviation using the SNR and maxampl
	ratio := math.Pow(10.0, float64(snr)/10.0)
	noiseSD = math.Sqrt(0.5 * float64(maxampl) * float64(maxampl) / ratio)
	file, err := os.Create(path.Join(dataDir, "sinsum.txt"))
	if err != nil {
		return fmt.Errorf("create %v error: %v", path.Join(dataDir, "sinsum.txt"), err)
	}
	defer file.Close()

	// Sum the sinsuoids and noise over the time interval
	delta := 1.0 / float64(samplerate)
	twoPi := 2.0 * math.Pi
	for t := 0.0; t < float64(samples)*delta; t += delta {
		sinesum := 0.0
		for _, sig := range signal {
			omega := twoPi * float64(sig.freq)
			sinesum += float64(sig.amp) * math.Sin(omega*t)
		}
		sinesum += noiseSD * rand.NormFloat64()
		fmt.Fprintf(file, "%v %v\n", t, sinesum)
	}

	return nil
}

// generateAM creates an amplitude modulation waveform
func generateAM(w http.ResponseWriter, r *http.Request, samples int, samplerate int) error {
	// get modulation percent
	temp := r.FormValue("percentmodulation")
	if len(temp) == 0 {
		return fmt.Errorf("missing Modulation Percent for AM")
	}
	percentModulation, err := strconv.Atoi(temp)
	if err != nil {
		return err
	}

	// get carrier frequency
	temp = r.FormValue("carrierfreq")
	if len(temp) == 0 {
		return fmt.Errorf("missing Carrier Frequency  for AM")
	}
	f1, err := strconv.Atoi(temp)
	if err != nil {
		return err
	}

	// get modulating frequency
	temp = r.FormValue("modulfreq")
	if len(temp) == 0 {
		return fmt.Errorf("missing Modulating Frequency for AM")
	}
	f2, err := strconv.Atoi(temp)
	if err != nil {
		return err
	}

	// alias check
	if f1-f2 < 0 {
		return fmt.Errorf("aliasing:  carrier frequency - modulating frequency < 0")
	}
	if f1+f2 > samplerate/2 {
		return fmt.Errorf("aliasing:  carrier frequency + modulating frequency > samplerate/2")
	}

	file, err := os.Create(path.Join(dataDir, "am.txt"))
	if err != nil {
		return fmt.Errorf("create %v error: %v", path.Join(dataDir, "am.txt"), err)
	}
	defer file.Close()

	// Create the AM waveform using sin(2*PI*f1*t)*(1 + A*cos(2*PI*f2*t))
	// where f1 is the carrier frequency, f2 is the modulating frequency,
	// and A is the percent modulation/100
	delta := 1.0 / float64(samplerate)
	const twoPi = 2.0 * math.Pi
	A := float64(percentModulation) / 100.0
	for t := 0.0; t < float64(samples)*delta; t += delta {
		y1 := math.Sin(twoPi * float64(f1) * t)
		y2 := 1.0 + A*math.Cos(twoPi*float64(f2)*t)
		fmt.Fprintf(file, "%v %v\n", t, y1*y2)
	}

	return nil
}

// generateFM creates a frequency modulation waveform, either sine baseband or linear (LFM)
func generateFM(w http.ResponseWriter, r *http.Request, samples int, samplerate int) error {

	generateLFM := func() error {
		// get bandwidth
		temp := r.FormValue("FMbandwidth")
		if len(temp) == 0 {
			return fmt.Errorf("missing Bandwidth  for FM linear")
		}
		bw, err := strconv.Atoi(temp)
		if err != nil {
			return err
		}

		// get center frequency
		temp = r.FormValue("FMfrequency")
		if len(temp) == 0 {
			return fmt.Errorf("missing Center Frequency for FM linear")
		}
		fc, err := strconv.Atoi(temp)
		if err != nil {
			return err
		}

		// check for aliasing
		if fc-bw/2 < 0 {
			return fmt.Errorf("aliasing:  fc-bw/2 < 0")
		}
		if fc+bw/2 > samplerate/2 {
			return fmt.Errorf("aliasing: fc+bw/2 > samplerate/2")
		}

		file, err := os.Create(path.Join(dataDir, "fm.txt"))
		if err != nil {
			return fmt.Errorf("create %v error: %v", path.Join(dataDir, "fmlinear.txt"), err)
		}
		defer file.Close()

		// Create the FM waveform using A*cos(2*PI*(fc-bw/2)*t + PI*(bw/tau)*t*t)
		// where fc is the center frequency, bw is the bandwidth, tau is the pulse
		// duration in seconds
		delta := 1.0 / float64(samplerate)
		tau := delta * float64(samples)
		const (
			pi    = math.Pi
			twoPi = 2.0 * pi
		)
		A := 1.0
		for t := 0.0; t < tau; t += delta {
			ang := twoPi*(float64(fc)-float64(bw)/2.0)*t + pi*(float64(bw)/tau)*t*t
			fmt.Fprintf(file, "%v %v\n", t, A*math.Cos(ang))
		}
		return nil
	}

	generateSineFM := func() error {
		// get center frequency
		temp := r.FormValue("FMfrequency")
		if len(temp) == 0 {
			return fmt.Errorf("missing Center Frequency for FM sine")
		}
		fc, err := strconv.Atoi(temp)
		if err != nil {
			return err
		}

		// get frequency deviation
		temp = r.FormValue("FMfreqdev")
		if len(temp) == 0 {
			return fmt.Errorf("missing frequency deviation for FM sine")
		}
		fd, err := strconv.Atoi(temp)
		if err != nil {
			return err
		}

		// get modulating frequency
		temp = r.FormValue("FMmodfreq")
		if len(temp) == 0 {
			return fmt.Errorf("missing modulating frequency for FM sine")
		}
		fm, err := strconv.Atoi(temp)
		if err != nil {
			return err
		}

		// check for aliasing
		if fc-fd < 0 {
			return fmt.Errorf("aliasing:  fc-bw/2 < 0")
		}
		if fc+fd > samplerate/2 {
			return fmt.Errorf("aliasing: fc+bw/2 > samplerate/2")
		}

		file, err := os.Create(path.Join(dataDir, "fm.txt"))
		if err != nil {
			return fmt.Errorf("create %v error: %v", path.Join(dataDir, "fmsine.txt"), err)
		}
		defer file.Close()

		// Create the FM waveform using A*cos(2*PI*fc*t + (fd/fm)*sin(2*PI*fm*t))
		// where fc is the center frequency, fd is the frequency deviation, fm is
		// the modulating frequency.
		delta := 1.0 / float64(samplerate)
		tau := delta * float64(samples)
		const (
			twoPi = 2.0 * math.Pi
		)
		A := 1.0
		for t := 0.0; t < tau; t += delta {
			y := A * math.Cos(twoPi*float64(fc)*t+
				float64(fd)/float64(fm)*math.Sin(twoPi*float64(fm)*t))
			fmt.Fprintf(file, "%v %v\n", t, y)
		}

		return nil
	}

	// Determine if Sine baseband or linear (LFM) and create it
	temp := r.FormValue("FMtype")
	switch temp {
	case "linear":
		return generateLFM()
	case "sine":
		return generateSineFM()
	}
	return nil
}

// processSignalData generates and saves to disk the desired signal
func processSignalData(w http.ResponseWriter, r *http.Request, samples int, samplerate int) error {

	// Determine which signal to generate
	switch r.FormValue("signaltype") {
	case "sumsinusoids":
		err := generateSinsuoids(w, r, samples, samplerate)
		if err != nil {
			return err
		}
	case "am":
		err := generateAM(w, r, samples, samplerate)
		if err != nil {
			return err
		}
	case "fm":
		err := generateFM(w, r, samples, samplerate)
		if err != nil {
			return err
		}
	}

	return nil
}

// handleGenerateData creates figure or signal data and saves to disk
func handleGenerateData(w http.ResponseWriter, r *http.Request) {
	var plot PlotT

	samplestxt := r.FormValue("samples")
	sampleratetxt := r.FormValue("samplerate")
	// choose time or frequency domain processing
	if len(samplestxt) > 0 && len(sampleratetxt) > 0 {

		samples, err := strconv.Atoi(samplestxt)
		if err != nil {
			plot.Status = fmt.Sprintf("Samples conversion to int error: %v", err.Error())
			fmt.Printf("Samples conversion to int error: %v", err.Error())
			// Write to HTTP using template and grid
			if err := generateTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		samplerate, err := strconv.Atoi(sampleratetxt)
		if err != nil {
			plot.Status = fmt.Sprintf("Sample rate conversion to int error: %v", err.Error())
			fmt.Printf("Sample rate conversion to int error: %v", err.Error())
			// Write to HTTP using template and grid
			if err := generateTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		// Get the category of data to generate
		switch r.FormValue("category") {
		case "figure":
			err := processFigureData(w, r, samples, samplerate)
			if err != nil {
				plot.Status = err.Error()
				if err := generateTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}
		case "signal":
			err := processSignalData(w, r, samples, samplerate)
			if err != nil {
				plot.Status = err.Error()
				if err := generateTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}
		}
		plot.Status = "Figures or Signals files written to the data folder"
		// Write to HTTP using template and grid
		if err := generateTmpl.Execute(w, plot); err != nil {
			log.Fatalf("Write to HTTP output using template with error: %v\n", err)
		}

	} else {
		plot.Status = "Enter number of samples and sample rate"
		// Write to HTTP using template and grid
		if err := generateTmpl.Execute(w, plot); err != nil {
			log.Fatalf("Write to HTTP output using template with error: %v\n", err)
		}
	}
}

// executive program
func main() {
	// Setup http server with handler for reading file time data, calculating frequency domain data,
	// and plotting time or frequency domain data
	http.HandleFunc(patternPlotData, handlePlotData)

	// Setup http server with handler for generating data
	http.HandleFunc(patternGenerateData, handleGenerateData)

	// Setup http server with handler for choosing data file and either time or frequency domain plot
	http.HandleFunc(patternPlotOptions, handlePlotOptions)
	fmt.Printf("Plot Data Server listening on %v.\n", addr)

	http.ListenAndServe(addr, nil)
}
