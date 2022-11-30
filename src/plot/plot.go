/*
Plot data in x-y coordinates from a file, one space-separated coordinate per line.
x1 y1
x2 y2
...
xn yn

Use CSS display grid to display a 100x100 grid of cells.
*/

package main

import (
	"bufio"
	"fmt"
	"log"
	"net/http"
	"os"
	"strconv"
	"strings"
	"text/template"
)

const (
	rows             = 100                                      // #rows in grid
	columns          = 100                                      // #columns in grid
	tmpl             = "../../src/plot/templates/plotdata.html" // html template relative address
	addr             = "127.0.0.1:8080"                         // http server listen address
	pattern          = "/plotdata"                              // http handler pattern for plotting data
	patternGenerator = "/generatedata"                          // http handler pattern for data generation
	xlabels          = 10                                       // # labels on x axis
)

type PlotT struct {
	Grid   []string // plotting grid
	Status string   // status of the plot
	Ymax   string
	Ymin   string
	Xlabel []string
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

	filename := r.FormValue("filename")
	if len(filename) > 0 {
		// Open file
		f, err := os.Open(filename)
		if err == nil {
			// Mark the data x-y coordinate online at the corresponding
			// grid row/column.
			input := bufio.NewScanner(f)

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
				} else if x < xmin {
					xmin = x
				}

				if y > ymax {
					ymax = y
				} else if y < ymin {
					ymin = y
				}
			}
			f.Close()
			f, err = os.Open(filename)
			if err == nil {
				input := bufio.NewScanner(f)

				// Calculate scale factors for x and y
				xscale = columns / (xmax - xmin)
				yscale = rows / (ymax - ymin)

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

					// This cell location (row,col) is on the line
					row := int(rows - (y-ymin)*yscale + .5)
					col := int((x-xmin)*xscale + .5)
					plot.Grid[row*rows+col] = "online"

				}
			} else {
				// Set plot status
				fmt.Printf("Error opening file %s: %v\n", filename, err)
				plot.Status = "Status: Error opening file."
			}

			// Set plot status
			plot.Status = fmt.Sprintf("Status: Data plotted from %s", filename)

		} else {
			// Set plot status
			fmt.Printf("Error opening file %s: %v\n", filename, err)
			plot.Status = "Status: Error opening file."
		}
	} else {
		plot.Status = "Status: Provide a data file to plot."
	}

	// Assign Ymax and Ymin
	plot.Ymax = fmt.Sprintf("%.2f", ymax)
	plot.Ymin = fmt.Sprintf("%.2f", ymin)

	// Construct x-axis labels
	incr := (xmax - xmin) / (xlabels - 1)
	x := xmin
	for i := range plot.Xlabel {
		plot.Xlabel[i] = fmt.Sprintf("%.2f", x)
		x += incr
	}

	// Write to HTTP output using template and grid
	if err := t.Execute(w, plot); err != nil {

		log.Fatalf("Write to HTTP output using template with grid error: %v\n", err)
	}

}

// handleGenerating creates x-y data and saves to disk files
func handleGenerating(w http.ResponseWriter, r *http.Request) {

}

// executive program
func main() {
	// Setup http server with handler for reading file data and plotting it
	http.HandleFunc(pattern, handlePlotting)
	// Setup http server with handler for generating data for testing
	http.HandleFunc(patternGenerator, handleGenerating)
	http.ListenAndServe(addr, nil)
}
