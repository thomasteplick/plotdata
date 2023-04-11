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
	"time"

	"github.com/mjibson/go-dsp/fft"
)

const (
	rows               = 300                                // #rows in grid
	columns            = 300                                // #columns in grid
	tmpltime           = "templates/plottimedata.html"      // html template address
	tmplfrequency      = "templates/plotfrequencydata.html" // html template address
	tmploptions        = "templates/plotoptions.html"       // html template address
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
	timeTmpl   *template.Template
	freqTmpl   *template.Template
	optionTmpl *template.Template
	winType    = []string{"Bartlett", "Welch", "Hamming", "Hanning"}
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

// init parses the html template files are done only once
func init() {
	timeTmpl = template.Must(template.ParseFiles(tmpltime))
	freqTmpl = template.Must(template.ParseFiles(tmplfrequency))
	optionTmpl = template.Must(template.ParseFiles(tmploptions))
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

	// Write to HTTP using template and grid``
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
		N = int(math.Exp2(math.Log2(float64(N) + .5)))

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
