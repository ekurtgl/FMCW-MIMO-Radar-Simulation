function [X,Y,Z] = ellipsoid2P(P1,P2,a,b,c,N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Construct an ellipsoid connecting two center points: P1 and P2.
% Semi-axis lengths are a, b, and c.
% N is the number of grid points for plotting the ellipsoid. 
% 
% By V.C. Chen
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cntr = (P1+P2)/2;  % ellipsoid center
Lc = norm(P2-P1);

% the axis defined by: P1+V*[0:Lc]
V = (P1-P2)/Lc;   %normalized cylinder's axis-vector;
U = rand(1,3);     %linear independent vector
U = V-U/(U*V');    %orthogonal vector to V
U = U/sqrt(U*U');  %orthonormal vector to V
W = cross(V,U);    %vector orthonormal to V and U
W = W/sqrt(W*W');  %orthonormal vector to V and U 

% generate the ellipsoid at (0,0,0)
[Xc,Yc,Zc] = ellipsoid(0,0,0,a,b,c,N);

A = kron(U',Xc);
B = kron(W',Yc);
C = kron(V',Zc);
TMP = A+B+C;
nt = size(TMP,2);

X = TMP(1:nt,:)+Cntr(1);
Y = TMP(nt+1:2*nt,:)+Cntr(2);
Z = TMP(2*nt+1:end,:)+Cntr(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
