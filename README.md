# plotdata
Create figures or signal waveforms and plot the data in the time or frequency domains with zoom and rotation features.  
 This is a web application that is written in Go that makes use of the standard library package http/template to generate the HTML.
 The http server is started with "go run plot.go" from the Command Prompt.  In a web browser Address bar, enter "http://127.0.0.1:8080/generatedata" or
 "http://127.0.0.1:8080/plotdata" to generate or plot data, respectively.  For generating data, various figures or signal waveforms can be created.  They
 will be saved to a disk file in the data folder.  The data can be plotted in the time or frequency domains.  For time domain plots, the data can be rotated
 or zoomed into to inspect in detail areas of interest.  In the frequency domain, periodograms are averaged with 50% overlapping segments.  The data is 
 windowed to reduce spectral leakage at a cost of reduced resolution because of the wider window main lobe.  The user can select the FFT size, number of 
 overlapped segments, log or linear plot, and window type.  The figures that can be created are a sine, circle, cardoid, sawtooth, and spiral.  The signals that 
 can be created are:
 <ul>
 <li>one to five sinsuoids in Gaussian noise at a specified signal-to-noise ratio (SNR)</li>
 <li>amplitude modulated (AM) waveform with desired center frequency and bandwidth</li>
 <li>frequency modulated (FM) waveform using linear frequency modulation (LFM) or baseband sine at desired center frequency, bandwidth, frequency deviation, or modulating frequency</li>
 </ul>
 
 The user can specify the number of samples and the sampling frequency used to create a figure or signal.

<h3>Generate Data Options Page</h3>

![image](https://user-images.githubusercontent.com/117768679/233181558-f3e2e0de-7356-4791-afa8-4c3dda775345.png)
<h3>Plot Data Options Page</h3>

![image](https://user-images.githubusercontent.com/117768679/233182057-4b6cba6b-9ed1-45c6-9a0d-a2aab2f9c8e6.png)

<h3>Spiral Figure in Time Domain</h3>

![image](https://user-images.githubusercontent.com/117768679/233192875-7d56a43c-07d9-42ae-b16e-0af15e8554fb.png)

<h3>Sine Figure in Time Domain</h3>

![image](https://user-images.githubusercontent.com/117768679/233193411-704f1c35-f0b0-4a28-b143-1af68c1d3559.png)

<h3>Five Sinsouoids in Gaussian Noise with 20 db SNR
 
 ![image](https://user-images.githubusercontent.com/117768679/233195762-c28a07c0-93ae-41cf-964d-40b4b8120e50.png)

 <h3>Five Sinsuoids in Gaussian Noise with 20 db SNR
  
  ![image](https://user-images.githubusercontent.com/117768679/233196675-4dcaeb7d-d5fc-430e-8f54-91faf8924fe8.png)

  <h3>Linear Frequency Modulation, 1000 Hz center frequency, 200 Hz bandwidth
   
   ![image](https://user-images.githubusercontent.com/117768679/233197696-cc9bab1b-9686-43c0-a905-f27bae158520.png)
   
   <h3>Linear Frequency Modulation, 1000 Hz center frequency, 400 Hz bandwidth
    
   ![image](https://user-images.githubusercontent.com/117768679/233197945-965151f3-adff-42d3-a678-92329be6ba2e.png)
 
