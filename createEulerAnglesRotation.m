function mat = createEulerAnglesRotation(phi, theta, psi)
%CREATEEULERANGLESROTATION Create a rotation matrix from 3 euler angles
%
%   ROT = createEulerAnglesRotation(PHI, THETA, PSI)
%   Create a rotation matrix from the 3 euler angles PHI THETA and PSI,
%   in radians, using the 'XYZ' convention. These angles correspond to the
%   "Roll-Pitch-Yaw" convention, also known as "TaitBryan angles".
%   PHI:    rotation angle around X-axis, in radians, corresponding to the
%       'Roll'. PHI is between -pi and +pi. 
%   THETA:  rotation angle around Y-axis, in radians, corresponding to the
%       'Pitch'. THETA is between -pi/2 and pi/2.
%   PSI:    rotation angle around Z-axis, in radians, corresponding to the
%       'Yaw'. PSI is between -pi and +pi.
%
%   The resulting rotation is equivalent to a rotation around X-axis by an
%   angle PHI, followed by a rotation around the Y-axis by an angle THETA,
%   followed by a rotation around the Z-axis by an angle PSI. 
%   That is:
%       ROT = Rz * Ry * Rx;
%
%   Example
%   [n e f] = createCube;
%   phi     = 20*pi/180;
%   theta   = 30*pi/180;
%   psi     = 10*pi/180;
%   rot = createEulerAnglesRotation(phi, theta, psi);
%   n2 = transformPoint3d(n, rot);
%   drawPolyhedron(n2, f);
%
%   See also
%   transforms3d, createRotationOx, createRotationOy, createRotationOz
%   rotation3dAxisAndAngle
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2010-07-22,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

%   HISTORY
%   2011-06-20 deprecate

warning('MatGeom:deprecation', ...
    'Deprecated function, use ''eulerAnglesToRotation3d'' instead');

% create individual rotation matrices
rotX = createRotationOx(phi);
rotY = createRotationOy(theta);
rotZ = createRotationOz(psi);

% concatenate matrices
mat = rotZ * rotY * rotX;