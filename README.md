# FMCW-MIMO-Radar-Simulation

## This repository is a simulation frequency modulated continuous wave (FMCW), multiple input multiple output (MIMO) radars.

`FMCW_simulation.m` is the main script which creates point targets at estimates their range, velocity and angle information.

TX (blue) - RX (red) antenna positions:

![fig1, align=center](https://user-images.githubusercontent.com/66868163/185649842-fd3723d2-e444-442b-b846-ab9e09415674.png)

Transmitted, received and downconverted signals:

![fig2](https://user-images.githubusercontent.com/66868163/185650072-9e99732b-bcda-4c72-8d6d-b8e04a132595.png)

Range-Doppler Map:

![fig3](https://user-images.githubusercontent.com/66868163/185650125-66192f23-ec3c-401e-aa9c-21f1ff87467e.png)

Cell-averaging constant-false-alarm-rate (CA-CFAR) detection result:

![fig4](https://user-images.githubusercontent.com/66868163/185650285-46981394-df12-4500-9448-015fc4679b5d.png)

Angle estimation spectrum result with Fast Fourier Transform (FFT):

![fig6](https://user-images.githubusercontent.com/66868163/185650616-58612191-d7fc-4ffa-ad4f-a7de2b9c9c1c.png)

Range-angle map with Fast Fourier Transform (FFT):

![fig5](https://user-images.githubusercontent.com/66868163/185650425-9f06c596-939d-43e3-920b-d4270235a3fc.png)

Angle estimation spectrum result with MUltiple SIgnal Classification (MUSIC):

![fig7](https://user-images.githubusercontent.com/66868163/185650774-a7b01bc5-43ca-470d-a083-bd459377db12.png)

Range-angle map with MUltiple SIgnal Classification (MUSIC):

![fig9](https://user-images.githubusercontent.com/66868163/185650845-ea30ccdf-e0f7-42bf-ac7d-5f3d648b4f93.png)

Point cloud generation from estimated 3-D coordinates:

![fig8](https://user-images.githubusercontent.com/66868163/185650990-6c3112cd-06f3-4d15-916b-34cf28eb0e5e.png)

----------

`FMCW_sim_v2.m` is the main script which reads 3-D coordinates of human skeleton joints captured by Kinect v2 device. For more information about the data capture and skeleton extraction using Kinect, check out the [Kinect repository](https://github.com/ekurtgl/Kinect). This script has not finished yet.

