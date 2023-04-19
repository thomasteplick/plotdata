# plotdata
Create figures or signal waveforms and plot the data in the time or frequency domains with zoom and rotation features.  
 This is a web application that is written in Go that makes use of the standard library package http/template to generate the HTML.
 The http server is started with "go run plot.go" from the Command Prompt.  In a web browser Address bar, enter "http://127.0.0.1:8080/generatedata" or
 "http://127.0.0.1:8080/plotdata" to generate or plot data, respectively.  For generating data, various figures or signal waveforms can be created.  They
 will be saved to a disk file in the data folder.  The data can be plotted in the time or frequency domains.  For time domain plots, the data can be rotated
 or zoomed into to inspect in detail areas of interest.
 
