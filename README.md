# plotdata
This program is a web application written in Go.  Build the Go source code or use the Run command to run the server in a Window's command shell .  In your web browser, 
enter http://127.0.0.1:8080/plotdata in the address box. Data files can be placed in the data/ folder which is in the same folder as plot.go.  Samples data files
are already located there.  The data in the text file are space separated numbers, two numbers per line.  The data can be zoomed into to inspect sections of interest.
The plots can also be rotated to any desired angle in degrees.

The display is a 300 x 300 cell grid, each cell is 2px by 2px.  Therefore there are 90,000 cells or 600 x 600 pixels = 360,000 pixels.  Extensive use is made of the
html/template package.

Entering http://127.0.0.1:8080/generatedata in the web browser will generate the text data files located in the data/ folder.

![image](https://user-images.githubusercontent.com/117768679/208189803-cb638b3a-8a10-4914-9041-131ece96e337.png)
![image](https://user-images.githubusercontent.com/117768679/208265667-91fad1c0-534b-4e1c-9d42-30c470c0f511.png)
![image](https://user-images.githubusercontent.com/117768679/208265736-50d0e597-6080-4455-ae9f-7ffb26b0358b.png)
![image](https://user-images.githubusercontent.com/117768679/208752632-c9f15bbd-ff1c-40ee-8cb8-807d82b03a27.png)
