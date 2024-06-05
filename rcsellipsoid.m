function rcs = rcsellipsoid(a,b,c,phi,theta) 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% This program calculates the backscattered RCS of an ellipsoid 
% with semi-axis lengths of a, b, and c. 
% The source code is based on Radar Systems Analysis and Design Using 
% MATLAB, By B. Mahafza, Chapman & Hall/CRC 2000. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
rcs = (pi*a^2*b^2*c^2)/(a^2*(sin(theta))^2*(cos(phi))^2+... 
    b^2*(sin(theta))^2*(sin(phi))^2+c^2*(cos(theta))^2)^2; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
