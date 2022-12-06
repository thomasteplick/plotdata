/*
Plot data in x-y coordinates from a file, one space-separated coordinate per line.
x1 y1
x2 y2
...
xn yn

The html/template package is used to generate the html sent to the client.
Use CSS display grid to display a 100x100 grid of cells.
Use CSS flexbox to display the labels on the x and y axes.
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
	rows             = 200                                      // #rows in grid
	columns          = 200                                      // #columns in grid
	tmpl             = "../../src/plot/templates/plotdata.html" // html template relative address
	addr             = "127.0.0.1:8080"                         // http server listen address
	pattern          = "/plotdata"                              // http handler pattern for plotting data
	patternGenerator = "/generatedata"                          // http handler pattern for data generation
	xlabels          = 11                                       // # labels on x axis
	ylabels          = 11                                       // # labels on y axis
	dataDir          = "data/"                                  // directory for the data files
)

type PlotT struct {
	Grid   []string // plotting grid
	Status string   // status of the plot
	Xlabel []string // x-axis labels
	Ylabel []string // y-axis labels
}

var (
	t *template.Template
)

// init parses the html template file done only once
func init() {
	t = template.Must(template.ParseFiles(tmpl))
}

// handlePlotting opens the file and plots the data
func handlePlotting(w http.ResponseWriter, r *http.Request) {
	// main data structure
	var (
		plot   PlotT
		xmax   float64
		xmin   float64
		ymax   float64
		ymin   float64
		xscale float64
		yscale float64
	)

	plot.Grid = make([]string, rows*columns)
	plot.Xlabel = make([]string, xlabels)
	plot.Ylabel = make([]string, ylabels)

	filename := r.FormValue("filename")
	if len(filename) > 0 {
		// Open file
		f, err := os.Open(filename)
		if err == nil {
			// Mark the data x-y coordinate online at the corresponding
			// grid row/column.
			input := bufio.NewScanner(f)
			xmax = -math.MaxFloat64
			xmin = math.MaxFloat64
			ymax = -math.MaxFloat64
			ymin = math.MaxFloat64

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

				if x > xmax {
					xmax = x
				}
				if x < xmin {
					xmin = x
				}

				if y > ymax {
					ymax = y
				}
				if y < ymin {
					ymin = y
				}
			}
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
						if (x1 < xmin || x1 > xmax) || (x2 < xmin || x2 > xmax) || (x1 >= x2) {
							plot.Status = "Zoom values are not in x range."
							fmt.Printf("Zoom error: start or end value not in x range.\n")
						} else if (y1 < ymin || y1 > ymax) || (y2 < ymin || y2 > ymax) || (y1 >= y2) {
							plot.Status = "Zoom values are not in y range."
							fmt.Printf("Zoom error: start or end value not in y range.\n")
						} else {
							// Valid Zoom endpoints, replace the previous min and max values
							xmin = x1
							xmax = x2
							ymin = y1
							ymax = y2
						}
					}
				}

				// Calculate scale factors for x and y
				xscale = (columns - 1) / (xmax - xmin)
				yscale = (rows - 1) / (ymax - ymin)

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

					// Check if inside the zoom values
					if x < xmin || x > xmax || y < ymin || y > ymax {
						continue
					}

					// This cell location (row,col) is on the line
					row := (rows - 1) - int((y-ymin)*yscale+.5)
					col := int((x-xmin)*xscale + .5)
					plot.Grid[row*columns+col] = "online"
				}
				// Set plot status if no errors
				if len(plot.Status) == 0 {
					plot.Status = fmt.Sprintf("Status: Data plotted from %s", filename)
				}

			} else {
				// Set plot status
				fmt.Printf("Error opening file %s: %v\n", filename, err)
				plot.Status = "Status: Error opening file."
			}

		} else {
			// Set plot status
			fmt.Printf("Error opening file %s: %v\n", filename, err)
			plot.Status = "Status: Error opening file."
		}
	} else {
		plot.Status = "Status: Provide a data file to plot."
	}

	// Construct x-axis labels
	incr := (xmax - xmin) / (xlabels - 1)
	x := xmin
	// First label is empty for alignment purposes
	for i := range plot.Xlabel {
		plot.Xlabel[i] = fmt.Sprintf("%.2f", x)
		x += incr
	}

	// Construct the y-axis labels
	incr = (ymax - ymin) / (ylabels - 1)
	y := ymin
	for i := range plot.Ylabel {
		plot.Ylabel[i] = fmt.Sprintf("%.2f", y)
		y += incr
	}

	// Write to HTTP using template and grid
	if err := t.Execute(w, plot); err != nil {
		log.Fatalf("Write to HTTP output using template with grid error: %v\n", err)
	}
}

// handleGenerating creates x-y data and saves to disk files
func handleGenerating(w http.ResponseWriter, r *http.Request) {
	const (
		samples         = 800
		cycles          = 4
		k       float64 = 2.0 * math.Pi / (samples / cycles)
	)

	// 1. sine wave
	f, err := os.Create(path.Join(dataDir, "sin.txt"))
	if err != nil {
		log.Fatalf("Create %v error: %v", path.Join(dataDir, "sin.txt"), err)
	}

	for i := 0; i < samples; i++ {
		fmt.Fprintf(f, "%v %v\n", i, math.Sin(k*float64(i)))
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
	// Setup http server with handler for reading file data and plotting it
	http.HandleFunc(pattern, handlePlotting)
	// Setup http server with handler for generating data for testing
	http.HandleFunc(patternGenerator, handleGenerating)
	http.ListenAndServe(addr, nil)
}
