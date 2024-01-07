# FMCW-MIMO-Radar-Simulation

## This repository is a simulation of frequency modulated continuous wave (FMCW), multiple input multiple output (MIMO) radars.

`FMCW_simulation.m` is the main script which creates point targets and estimates their range, velocity and angle information.

TX (blue) - RX (red) antenna positions:

<p align="center">
  <img width="460" height="300" src="https://user-images.githubusercontent.com/66868163/185649842-fd3723d2-e444-442b-b846-ab9e09415674.png">
</p>

Transmitted, received and downconverted signals:

<p align="center">
  <img width="460" height="300" src="https://user-images.githubusercontent.com/66868163/185650072-9e99732b-bcda-4c72-8d6d-b8e04a132595.png">
</p>

Range-Doppler Map:

<p align="center">
  <img width="460" height="300" src="https://user-images.githubusercontent.com/66868163/185650125-66192f23-ec3c-401e-aa9c-21f1ff87467e.png">
</p>

Cell-averaging constant-false-alarm-rate (CA-CFAR) detection result:

<p align="center">
  <img width="460" height="300" src="https://user-images.githubusercontent.com/66868163/185650285-46981394-df12-4500-9448-015fc4679b5d.png">
</p>

Angle estimation spectrum result with Fast Fourier Transform (FFT):

<p align="center">
  <img width="460" height="300" src="https://user-images.githubusercontent.com/66868163/185650616-58612191-d7fc-4ffa-ad4f-a7de2b9c9c1c.png">
</p>

Range-angle map with Fast Fourier Transform (FFT):

<p align="center">
  <img width="460" height="300" src="https://user-images.githubusercontent.com/66868163/185650425-9f06c596-939d-43e3-920b-d4270235a3fc.png">
</p>

Angle estimation spectrum result with MUltiple SIgnal Classification (MUSIC):

<p align="center">
  <img width="460" height="300" src="https://user-images.githubusercontent.com/66868163/185650774-a7b01bc5-43ca-470d-a083-bd459377db12.png">
</p>

Range-angle map with MUltiple SIgnal Classification (MUSIC):

<p align="center">
  <img width="460" height="300" src="https://user-images.githubusercontent.com/66868163/185650845-ea30ccdf-e0f7-42bf-ac7d-5f3d648b4f93.png">
</p>

Point cloud generation from estimated 3-D coordinates:

<p align="center">
  <img width="460" height="300" src="https://user-images.githubusercontent.com/66868163/185650990-6c3112cd-06f3-4d15-916b-34cf28eb0e5e.png">
</p>

----------

`FMCW_sim_v3.m` is the main script which reads 3-D coordinates of human skeleton joints captured by Kinect v2 device, and extracts the raw radar data cube (RDC) and plays range-Doppler maps, and outputs the micro-Doppler spectrogram. For more information about the data capture and skeleton extraction using Kinect, check out the [Kinect repository](https://github.com/ekurtgl/Kinect). `This script is still under development!`

Range-Doppler Map:

<p align="center">
  <img width="460" height="300" src="https://user-images.githubusercontent.com/66868163/185657077-0aa8f056-e375-4ff4-97c9-b69f897ba2c0.gif">
</p>

Micro-Doppler Spectrogram:

<p align="center">
  <img width="460" height="300" src="https://user-images.githubusercontent.com/66868163/186519037-11e4f640-7861-4522-a94c-ad7e95e813e0.png">
</p>

# Boulic Model
`FMCW_sim_v4.m` simulates the Boulic Human Walking Model:

<p align="center">
  <img width="460" height="300" src="https://github.com/ekurtgl/FMCW-MIMO-Radar-Simulation/assets/66868163/18f9ee8d-2f10-4377-b98b-9dc94250912b">
</p>

<p align="center">
  <img width="460" height="300" src="https://github.com/ekurtgl/FMCW-MIMO-Radar-Simulation/assets/66868163/ff820010-2fc6-4610-9d8f-8594c3c3c40b">
</p>

# Range-Elevation Map for Real VICON MOCAP Data

<p align="center">
  <img width="460" height="300" src="https://github.com/ekurtgl/FMCW-MIMO-Radar-Simulation/assets/66868163/4f74d752-82a6-43a3-8d81-5091b55efe27">
</p>
