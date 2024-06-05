function [X Y Z] = cylinder2P(P1,P2,R,N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Construct a cylinder connecting two center points: P1 and P2.
% R = [R1; R2], where R1 and R2 are radiuses of the cylinder's two ends.
%               If R1 = R2, R can be a variable.
% N is the number of grid points for plotting the cylinder. 
% 
% By V.C. Chen
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(R)==1
    R = [R;R];
end
        
% Calculating the length of the cylinder
Lc = norm(P2-P1);

% angle array
theta = linspace(0,2*pi,N)';

X = zeros(length(R), N); 
Y = zeros(length(R), N);
Z = zeros(length(R), N);

% cylinder axis definedd by: P1+V*[0:Lc]
V = (P2-P1)/Lc;   %normalized cylinder's axis-vector;
U = rand(1,3);    %linear independent vector
U = V-U/(U*V');   %orthogonal vector to V
U = U/sqrt(U*U'); %orthonormal vector to V
W = cross(V,U);   %vector orthonormal to V and U
W = W/sqrt(W*W'); %orthonormal vector to V and U 

P1x = P1(1);
P1y = P1(2);
P1z = P1(3);
P2x = P2(1);
P2y = P2(2);
P2z = P2(3);

Vx = V(1);
Vy = V(2);
Vz = V(3);

Ux = U(1);
Uy = U(2);
Uz = U(3);

Wx = W(1);
Wy = W(2);
Wz = W(3);

Nrad = linspace(0,1,2);

for k=1:length(R)
  r = Nrad(k);
  X(k, :) = P1x + (P2x-P1x)*r + R(k)*sin(theta)*Ux-R(k)*cos(theta)*Wx; 
  Y(k, :) = P1y + (P2y-P1y)*r + R(k)*sin(theta)*Uy-R(k)*cos(theta)*Wy; 
  Z(k, :) = P1z + (P2z-P1z)*r + R(k)*sin(theta)*Uz-R(k)*cos(theta)*Wz;  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
