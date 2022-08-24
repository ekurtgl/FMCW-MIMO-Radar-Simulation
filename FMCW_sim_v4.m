clear; clc; close all

load('20220817_191142_00.mat');
%% Radar parameters
c = physconst('LightSpeed'); %speed of light
BW = 4e9; %bandwidth
start_freq = 77e9;
end_freq = start_freq + BW;
fc = (start_freq + end_freq)/2; % carrier frequency
numADC = 256; % # of adc samples
% numChirps = 256; % # of chirps per frame
% numCPI = 10;
NPpF = 256;
powerTx = 30; % dB
frameDuration = 20e-3; % 40e-3
T = frameDuration/NPpF; % PRI
PRF = 1/T;
F = numADC/T; % sampling frequency
dt = 1/F; % sampling interval
slope = BW/T;
lambda = c/fc;
% N = numChirps*numADC*numCPI; % total # of adc samples
% t = linspace(0,T*numChirps*numCPI,N); % time axis, one frame
% t= 0:dt:dt*numADC*numChirps-dt;
t_onePulse = 0:dt:dt*numADC-dt;
% t_onePulse = dt:dt:dt*numADC;
% t_onePulse = linspace(0, T, numADC);
numTX = 1;
numRX = 1;
Vmax = lambda/(T*4); % Max Unamb velocity m/s
DFmax = 1/2*PRF; % = Vmax/(c/fc/2); % Max Unamb Dopp Freq
dR = c/(2*BW); % range resol
Rmax = F*c/(2*slope); % TI's MIMO Radar doc
Rmax2 = c/2/PRF; % lecture 2.3
% dV = lambda/(2*numChirps*T); % velocity resol, lambda/(2*framePeriod)
d_rx = lambda/2; % dist. between rxs
d_tx = 4*d_rx; % dist. between txs

% N_Dopp = numChirps; % length of doppler FFT
N_range = numADC; % length of range FFT
N_azimuth = numTX*numRX;
R = 0:dR:Rmax-dR; % range axis
% V = linspace(-Vmax, Vmax, numChirps); % Velocity axis
ang_ax = -90:90; % angle axis

%% Antennas
tx_loc = cell(1,numTX);
for i = 1:numTX
   tx_loc{i} = [(i-1)*d_tx 0 0];
%    scatter3(tx_loc{i}(1),tx_loc{i}(2),tx_loc{i}(3),'b','filled')
%    hold on
end

rx_loc = cell(1,numRX);
for i = 1:numRX
   rx_loc{i} = [tx_loc{numTX}(1)+d_tx+(i-1)*d_rx 0 0];
%    scatter3(rx_loc{i}(1),rx_loc{i}(2),rx_loc{i}(3),'r','filled')
end

%% Targets
fps_skel = 30;
num_tar = size(skel_hist,1);
durationx = size(skel_hist,3)/fps_skel;
numChirps = floor(durationx*NPpF*(1/frameDuration));
tar_loc = zeros(num_tar, size(skel_hist,2), numChirps);
vel_hist = zeros(num_tar, size(skel_hist,2));

for t = 1:num_tar
    for i = 1:size(skel_hist,2)
        tar_loc(t,i,:) = spline(1:size(skel_hist,3), squeeze(skel_hist(t,i,:)), linspace(1,size(skel_hist,3), numChirps));      
%         tar_loc(t,i,:) = makima(1:size(skel_hist,3), squeeze(skel_hist(t,i,:)), linspace(1,size(skel_hist,3), numChirps));      
        vel_hist(t,i) = (max(tar_loc(t,i,:)) - min(tar_loc(t,i,:))) * sqrt(3) / durationx;
    end
end
% tar_loc = tar_loc+2;
% tt = squeeze(tar_loc(1,3,:));
v_avg = mean(vel_hist,2);
% v_avg = (max(tar_loc(1,2,:)) - min(tar_loc(1,2,:))) * sqrt(3) / durationx;

% figure
% hold on
% t = 2;
% for i = 1:numADC*10:numChirps*numADC
%    scatter3(tar_loc(t,1,i), tar_loc(t,2,i), tar_loc(t,3,i))
% end

%% Elipsoid model

% figure('visible','off')
% cnt = 1;
% for n = 1:ceil(numChirps/size(skel_hist,3)):numChirps
%     disp(['Frame ' int2str(n) '/' int2str(numChirps)]);
%     create_body_model(squeeze(tar_loc(:,:,n)).');
%     title(['Elipsoid model, Chirp: ' int2str(n) '/' int2str(numChirps)]);
%     F2(cnt) = getframe(gcf); % gcf returns the current figure handle
%     cnt = cnt + 1;
% end
% 
% writerObj = VideoWriter('elipsoid.avi');
% writerObj.FrameRate = fps_skel;
% open(writerObj);
% 
% for i=1:length(F2)
%         frame = F2(i) ;
%         writeVideo(writerObj, frame);
% end
% close(writerObj);

% common variables
r_dist = zeros(3,numChirps);
distances = zeros(1,numChirps);
aspct = zeros(3,numChirps);
ThetaAngle = zeros(1,numChirps);
PhiAngle = zeros(1,numChirps);
rcs = zeros(1,numChirps);
amp = zeros(num_tar,numChirps);

% radar returns from the head
target_id = [3 4]; % head
head = zeros(3,numChirps);
neck = zeros(3,numChirps);
disp('Processing head ...')
for k = 1:numChirps
    
    head(:,k) = tar_loc(4,:,k);
    neck(:,k) = tar_loc(3,:,k);
    headlen = sqrt(sum((head(:,k)-neck(:,k)).^2));
    
    r_dist(:,k) = abs(head(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector of the head
    aspct(:,k) = head(:,k)-neck(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-head(1,k); tx_loc{1}(2)-head(2,k);...
        tx_loc{1}(3)-head(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-head(2,k))./...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.1; % ellipsoid parameter
    b = 0.1;
    c = headlen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((head(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
    
end

% radar returns from the torso
target_id = [1 2 21]; % spine base, spine mid
neck = zeros(3,numChirps);
base = zeros(3,numChirps);
torso = zeros(3,numChirps);
disp('Processing torso ...')
for k = 1:numChirps
    
    base(:,k) = tar_loc(1,:,k);
    neck(:,k) = tar_loc(3,:,k);
    torsolen = sqrt(sum((neck(:,k)-base(:,k)).^2));
    
    torso(:,k) = (neck(:,k)+base(:,k))/2;
    r_dist(:,k) = abs(torso(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector 
    aspct(:,k) = neck(:,k)-base(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-torso(1,k); tx_loc{1}(2)-torso(2,k);...
        tx_loc{1}(3)-torso(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-torso(2,k))./...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.15; % ellipsoid parameter
    b = 0.15;
    c = torsolen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((torso(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

% radar returns from the left shoulder
target_id = 5; 
lshoulder = zeros(3,numChirps);
disp('Processing left shoulder ...')
for k = 1:numChirps
    
    lshoulder(:,k) = tar_loc(5,:,k);
    neck(:,k) = tar_loc(3,:,k);
    shoulderlen = sqrt(sum((neck(:,k)-lshoulder(:,k)).^2));
    
    r_dist(:,k) = abs(lshoulder(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector 
    aspct(:,k) = lshoulder(:,k)-neck(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-lshoulder(1,k); tx_loc{1}(2)-lshoulder(2,k);...
        tx_loc{1}(3)-lshoulder(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-lshoulder(2,k))./...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.06;
    b = 0.06;
    c = shoulderlen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((lshoulder(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

% radar returns from the right shoulder
target_id = 9; 
rshoulder = zeros(3,numChirps);
disp('Processing right shoulder ...')
for k = 1:numChirps
    
    rshoulder(:,k) = tar_loc(9,:,k);
    neck(:,k) = tar_loc(3,:,k);
    shoulderlen = sqrt(sum((neck(:,k)-rshoulder(:,k)).^2));
    
    r_dist(:,k) = abs(rshoulder(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector 
    aspct(:,k) = rshoulder(:,k)-neck(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-rshoulder(1,k); tx_loc{1}(2)-rshoulder(2,k);...
        tx_loc{1}(3)-rshoulder(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-rshoulder(2,k))./...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.06;
    b = 0.06;
    c = shoulderlen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((rshoulder(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

% radar returns from the left upper arm
target_id = 6; 
lelbow = zeros(3,numChirps);
lupperarm = zeros(3,numChirps);
disp('Processing left upper arm ...')
for k = 1:numChirps
    
    lshoulder(:,k) = tar_loc(5,:,k);
    lelbow(:,k) = tar_loc(6,:,k);
    upperarmlen = sqrt(sum((lelbow(:,k)-lshoulder(:,k)).^2));
    
    lupperarm(:,k) = (lshoulder(:,k)+lelbow(:,k))/2;
    r_dist(:,k) = abs(lupperarm(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector
    aspct(:,k) = lshoulder(:,k)-lelbow(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-lupperarm(1,k); tx_loc{1}(2)-lupperarm(2,k);...
        tx_loc{1}(3)-lupperarm(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-lupperarm(2,k))./...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.06;
    b = 0.06;
    c = upperarmlen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((lupperarm(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

% radar returns from the right upper arm
target_id = 10; 
relbow = zeros(3,numChirps);
rupperarm = zeros(3,numChirps);
disp('Processing right upper arm ...')
for k = 1:numChirps
    
    rshoulder(:,k) = tar_loc(9,:,k);
    relbow(:,k) = tar_loc(10,:,k);
    upperarmlen = sqrt(sum((relbow(:,k)-rshoulder(:,k)).^2));
    
    rupperarm(:,k) = (rshoulder(:,k)+relbow(:,k))/2;
    r_dist(:,k) = abs(rupperarm(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector
    aspct(:,k) = rshoulder(:,k)-relbow(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-rupperarm(1,k); tx_loc{1}(2)-rupperarm(2,k);...
        tx_loc{1}(3)-rupperarm(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-rupperarm(2,k))./...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.06;
    b = 0.06;
    c = upperarmlen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((rupperarm(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

% radar returns from the left lower arm
target_id = 8; 
lhand = zeros(3,numChirps);
disp('Processing left lower arm ...')
for k = 1:numChirps
    
    lhand(:,k) = tar_loc(8,:,k);
    lelbow(:,k) = tar_loc(6,:,k);
    lowerarmlen = sqrt(sum((lelbow(:,k)-lhand(:,k)).^2));
    
    r_dist(:,k) = abs(lhand(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector
    aspct(:,k) = lelbow(:,k)-lhand(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-lhand(1,k); tx_loc{1}(2)-lhand(2,k);...
        tx_loc{1}(3)-lhand(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-lhand(2,k))./ ...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.05;
    b = 0.05;
    c = lowerarmlen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((lhand(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

% radar returns from the right lower arm
target_id = 12; 
rhand = zeros(3,numChirps);
disp('Processing right lower arm ...')
for k = 1:numChirps
    
    rhand(:,k) = tar_loc(12,:,k);
    relbow(:,k) = tar_loc(10,:,k);
    lowerarmlen = sqrt(sum((relbow(:,k)-rhand(:,k)).^2));
    
    r_dist(:,k) = abs(rhand(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector
    aspct(:,k) = relbow(:,k)-rhand(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-rhand(1,k); tx_loc{1}(2)-rhand(2,k);...
        tx_loc{1}(3)-rhand(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-rhand(2,k))./ ...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.05;
    b = 0.05;
    c = lowerarmlen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((rhand(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

% radar returns from the left hip
target_id = 13; 
lhip = zeros(3,numChirps);
disp('Processing left hip ...')
for k = 1:numChirps
    
    lhip(:,k) = tar_loc(13,:,k);
    base(:,k) = tar_loc(1,:,k);
    hiplen = sqrt(sum((lhip(:,k)-base(:,k)).^2));
    
    r_dist(:,k) = abs(lhip(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector
    aspct(:,k) = lhip(:,k)-base(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-lhip(1,k); tx_loc{1}(2)-lhip(2,k);...
        tx_loc{1}(3)-lhip(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-lhip(2,k))./ ...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.07;
    b = 0.07;
    c = hiplen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((lhip(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

% radar returns from the right hip
target_id = 17; 
rhip = zeros(3,numChirps);
disp('Processing right hip ...')
for k = 1:numChirps
    
    rhip(:,k) = tar_loc(17,:,k);
    base(:,k) = tar_loc(1,:,k);
    hiplen = sqrt(sum((rhip(:,k)-base(:,k)).^2));
    
    r_dist(:,k) = abs(rhip(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector
    aspct(:,k) = rhip(:,k)-base(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-rhip(1,k); tx_loc{1}(2)-rhip(2,k);...
        tx_loc{1}(3)-rhip(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-rhip(2,k))./ ...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.07;
    b = 0.07;
    c = hiplen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((rhip(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

% radar returns from the left upper leg
target_id = 14; 
lknee = zeros(3,numChirps);
lupperleg = zeros(3,numChirps);
disp('Processing left upper leg ...')
for k = 1:numChirps
    
    lhip(:,k) = tar_loc(13,:,k);
    lknee(:,k) = tar_loc(14,:,k);
    upperleglen = sqrt(sum((lhip(:,k)-lknee(:,k)).^2));
    
    lupperleg(:,k) = (lknee(:,k)+lhip(:,k))/2;
    r_dist(:,k) = abs(lupperleg(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector
    aspct(:,k) = lknee(:,k)-lhip(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-lupperleg(1,k); tx_loc{1}(2)-lupperleg(2,k);...
        tx_loc{1}(3)-lupperleg(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-lupperleg(2,k))./ ...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.07;
    b = 0.07;
    c = upperleglen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((lupperleg(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

% radar returns from the right upper leg
target_id = 18; 
rknee = zeros(3,numChirps);
rupperleg = zeros(3,numChirps);
disp('Processing right upper leg ...')
for k = 1:numChirps
    
    rhip(:,k) = tar_loc(17,:,k);
    rknee(:,k) = tar_loc(18,:,k);
    upperleglen = sqrt(sum((rhip(:,k)-rknee(:,k)).^2));
    
    rupperleg(:,k) = (rknee(:,k)+rhip(:,k))/2;
    r_dist(:,k) = abs(rupperleg(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector
    aspct(:,k) = rknee(:,k)-rhip(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-rupperleg(1,k); tx_loc{1}(2)-rupperleg(2,k);...
        tx_loc{1}(3)-rupperleg(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-rupperleg(2,k))./ ...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.07;
    b = 0.07;
    c = upperleglen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((rupperleg(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

% radar returns from the left lower leg
target_id = 15; 
lankle = zeros(3,numChirps);
llowerleg = zeros(3,numChirps);
disp('Processing left lower leg ...')
for k = 1:numChirps
    
    lankle(:,k) = tar_loc(15,:,k);
    lknee(:,k) = tar_loc(14,:,k);
    lowerleglen = sqrt(sum((lankle(:,k)-lknee(:,k)).^2));
    
    llowerleg(:,k) = (lankle(:,k)+lknee(:,k))/2;
    r_dist(:,k) = abs(llowerleg(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector
    aspct(:,k) = lankle(:,k)-lknee(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-llowerleg(1,k); tx_loc{1}(2)-llowerleg(2,k);...
        tx_loc{1}(3)-llowerleg(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-llowerleg(2,k))./ ...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.06;
    b = 0.06;
    c = lowerleglen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((llowerleg(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

% radar returns from the right lower leg
target_id = 19; 
rankle = zeros(3,numChirps);
rlowerleg = zeros(3,numChirps);
disp('Processing right lower leg ...')
for k = 1:numChirps
    
    rankle(:,k) = tar_loc(19,:,k);
    rknee(:,k) = tar_loc(18,:,k);
    lowerleglen = sqrt(sum((rankle(:,k)-rknee(:,k)).^2));
    
    rlowerleg(:,k) = (rankle(:,k)+rknee(:,k))/2;
    r_dist(:,k) = abs(rlowerleg(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector
    aspct(:,k) = rankle(:,k)-rknee(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-rlowerleg(1,k); tx_loc{1}(2)-rlowerleg(2,k);...
        tx_loc{1}(3)-rlowerleg(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-rlowerleg(2,k))./ ...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.06;
    b = 0.06;
    c = lowerleglen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((rlowerleg(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

% radar returns from the left foot
target_id = 16; 
ltoe = zeros(3,numChirps);
disp('Processing right left foot ...')
for k = 1:numChirps
    
    lankle(:,k) = tar_loc(15,:,k);
    ltoe(:,k) = tar_loc(16,:,k);
    footlen = sqrt(sum((lankle(:,k)-ltoe(:,k)).^2));
    
    r_dist(:,k) = abs(ltoe(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector
    aspct(:,k) = lankle(:,k)-ltoe(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-ltoe(1,k); tx_loc{1}(2)-ltoe(2,k);...
        tx_loc{1}(3)-ltoe(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-ltoe(2,k))./ ...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.05;
    b = 0.05;
    c = footlen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((ltoe(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

% radar returns from the right foot
target_id = 20; 
rtoe = zeros(3,numChirps);
disp('Processing right right foot ...')
for k = 1:numChirps
    
    rankle(:,k) = tar_loc(19,:,k);
    rtoe(:,k) = tar_loc(20,:,k);
    footlen = sqrt(sum((rankle(:,k)-rtoe(:,k)).^2));
    
    r_dist(:,k) = abs(rtoe(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector
    aspct(:,k) = rankle(:,k)-rtoe(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-rtoe(1,k); tx_loc{1}(2)-rtoe(2,k);...
        tx_loc{1}(3)-rtoe(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-rtoe(2,k))./ ...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.05;
    b = 0.05;
    c = footlen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((rtoe(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

% radar returns from the left hand
target_id = [7 22 23]; 
lhandtip = zeros(3,numChirps);
lwrist = zeros(3,numChirps);
disp('Processing left hand ...')
for k = 1:numChirps
    
    lhandtip(:,k) = tar_loc(22,:,k);
    lwrist(:,k) = tar_loc(7,:,k);
    handlen = sqrt(sum((lhandtip(:,k)-lwrist(:,k)).^2));
    
    r_dist(:,k) = abs(lhandtip(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector
    aspct(:,k) = lwrist(:,k)-lhandtip(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-lhandtip(1,k); tx_loc{1}(2)-lhandtip(2,k);...
        tx_loc{1}(3)-lhandtip(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-lhandtip(2,k))./ ...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.06;
    b = 0.02;
    c = handlen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((lhandtip(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

% radar returns from the right hand
target_id = [11 24 25]; 
rhandtip = zeros(3,numChirps);
rwrist = zeros(3,numChirps);
disp('Processing right hand ...')
for k = 1:numChirps
    
    rhandtip(:,k) = tar_loc(24,:,k);
    rwrist(:,k) = tar_loc(11,:,k);
    handlen = sqrt(sum((rhandtip(:,k)-rwrist(:,k)).^2));
    
    r_dist(:,k) = abs(rhandtip(:,k)-tx_loc{1}.');
    distances(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    % aspect vector
    aspct(:,k) = rwrist(:,k)-rhandtip(:,k);
    % calculate theta angle
    A = [tx_loc{1}(1)-rhandtip(1,k); tx_loc{1}(2)-rhandtip(2,k);...
        tx_loc{1}(3)-rhandtip(3,k)];
    B = [aspct(1,k); aspct(2,k); aspct(3,k)];
    A_dot_B = dot(A,B,1);
    A_sum_sqrt = sqrt(sum(A.*A,1));
    B_sum_sqrt = sqrt(sum(B.*B,1));
    ThetaAngle(k) = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
    PhiAngle(k) = asin((tx_loc{1}(2)-rhandtip(2,k))./ ...
        sqrt(r_dist(1,k).^2+r_dist(2,k).^2));
    
    a = 0.06;
    b = 0.02;
    c = handlen/2;
    rcs(k) = rcsellipsoid(a,b,c,PhiAngle(k),ThetaAngle(k));
    
    distTx = distances(k);
    distRx = sqrt(sum((rhandtip(:,k)-rx_loc{1}').^2));
    amp(target_id, k) = sqrt(rcs(k)*powerTx) / (distTx * distRx);
end

%% TX
c = physconst('LightSpeed'); %speed of light
delays_targets = cell(numTX,numRX,num_tar);

for t = 1:num_tar
    for i = 1:numTX
        for j = 1:numRX
            delays_targets{i,j,t} = (vecnorm(squeeze(tar_loc(t,:,:))-repmat(rx_loc{j},size(tar_loc,3),1).',2,1)+vecnorm(squeeze(tar_loc(t,:,:))-repmat(tx_loc{i},size(tar_loc,3),1).',2,1))/c;
%             delays_targets{i,j,t} = delays_targets{i,j,t}(1:1:end);
        end
    end
end
dd=delays_targets{3};
numChirps = length(dd);
%% Complex signal
phase = @(tx,fx) 2*pi*(fx.*tx+slope/2*tx.^2); % transmitted
phase2 = @(tx,fx,r,v) 2*pi*(2*fx*r/c+tx.*(2*fx*v/c + 2*slope*r/c)); % downconverted

phase_t = phase(t_onePulse,fc);

mixed = zeros(numTX,numRX,numChirps*numADC);
% new_mixed = zeros(numTX, numRX, numADC, numChirps);
for i = 1:numTX
    for j = 1:numRX
        disp(['Processing Channel: ' num2str(j) '/' num2str(numRX)]);
        for t = 1:num_tar
            disp([int2str(t) '/' int2str(num_tar)]);
            
            for k = 1:numChirps
                phase_tar = phase(t_onePulse-delays_targets{i,j,t}(k),fc); % received
                signal_tar((k-1)*numADC+1:k*numADC) = exp(1j*(phase_t - phase_tar));
%                 signal_tar = exp(1j*2*pi*...
%                     (pdist([squeeze(tar_loc(t,:,k));rx_loc{j}])+pdist([squeeze(tar_loc(t,:,k));tx_loc{i}]))/lambda);
%                 new_mixed(i,j,:,k) = squeeze(new_mixed(i,j,:,k)) + signal_tar;
            end
%             mixed(i,j,:) = squeeze(mixed(i,j,:)) + signal_tar.';
            mixed(i,j,:) = squeeze(mixed(i,j,:)) + repelem(amp(t,:).', numADC, 1) .* signal_tar.';
%             break
        end
    end
end

% figure
% subplot(3,1,1)
% p1 = plot(t, real(signal_t));
% title('TX')
% xlim([0 0.1e-4])
% xlabel('Time (sec)');
% ylabel('Amplitude');
% subplot(3,1,2)
% p2 = plot(t, real(signal_1));
% title('RX')
% xlim([0 0.1e-4])
% xlabel('Time (sec)');
% ylabel('Amplitude');
% subplot(3,1,3)
% p3 = plot(t,real(mixed{i,j}));
% title('Mixed')
% xlim([0 0.1e-4])
% xlabel('Time (sec)');
% ylabel('Amplitude');

%% Post processing - 2-D FFT

% RDC = reshape(cat(3,mixed{:}),numADC,numChirps*numCPI,numRX*numTX); % radar data cube
RDC = reshape(mixed,numADC,numChirps,numRX*numTX);
numCPI = floor(numChirps/NPpF);
RDMs = zeros(numADC,NPpF,numTX*numRX,numCPI);
for i = 1:numCPI
    RD_frame = RDC(:,(i-1)*NPpF+1:i*NPpF,:);
    RDMs(:,:,:,i) = fftshift(fft2(RD_frame,[],[]),2);
end

% figure
% colormap(jet(256))
% for f = 1:numCPI
%     imagesc([-Vmax Vmax], [0 Rmax], 20*log10(abs(RDMs(:,:,1,f))/max(max(abs(RDMs(:,:,1,f))))));
%     clim = get(gca,'clim');
%     caxis([clim(1)/2 0])
% %     title(['Frame ' int2str(f) '/' int2str(numCPI)]);
%     xlabel('Velocity (m/s)');
%     ylabel('Range (m)');
%     title(['Range-Doppler Map, Frame: ' int2str(f) '/' int2str(numCPI)]);
%     drawnow;
%     F2(f) = getframe(gcf); % gcf returns the current figure handle
%     pause(frameDuration)
% end
% 
% writerObj = VideoWriter('test.avi');
% writerObj.FrameRate = floor(1/frameDuration);
% open(writerObj);

% for i=1:length(F2)
%         frame = F2(i) ;
%         writeVideo(writerObj, frame);
% end
% close(writerObj);
        
%% Micro-Doppler spectrogram

rBin = 1:256;
nfft = 2^12;window = 256;noverlap = 200;shift = window - noverlap;
sx = myspecgramnew(sum(fft(RDC(rBin,:,:))),window,nfft,shift); % mti filter and IQ correction
% sx = myspecgramnew(sum(new_mixed(rBin,:)),window,nfft,shift); % mti filter and IQ correction
sx2 = abs(flipud(fftshift(sx,1)));
timeAxis = (1:numCPI)*frameDuration; % Time
freqAxis = linspace(-PRF/2,PRF/2,nfft); % Frequency Axis
fig = figure('visible','on');
colormap(jet(256));
% set(gca,'units','normalized','outerposition',[0,0,1,1]);
doppSignMTI = imagesc(timeAxis,[-PRF/2 PRF/2],20*log10(abs(sx2/max(sx2(:)))));
%     axis xy
%     set(gca,'FontSize',10)
title('micro-Doppler Spectrogram');
%     title(fOut(end-22:end-4))
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
caxis([-45 0]) % 40
set(gca, 'YDir','normal')
set(gcf,'color','w');
%     colorbar;
% axis([0 timeAxis(end) -prf/6 prf/6])
%     saveas(fig,[fOut(1:end-4) '.fig']);
% set(gca,'xtick',[],'ytick',[])
