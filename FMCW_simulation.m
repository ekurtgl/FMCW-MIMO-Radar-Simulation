clear; clc; close all

%% Radar parameters
c = physconst('LightSpeed'); %speed of light
BW = 150e6; %bandwidth
fc = 77e9; % carrier frequency
numADC = 256; % # of adc samples
numChirps = 256; % # of chirps per frame
numCPI = 10;
T = 10e-6; % PRI 
PRF = 1/T;
F = numADC/T; % sampling frequency
dt = 1/F; % sampling interval
slope = BW/T;
lambda = c/fc;
N = numChirps*numADC*numCPI; % total # of adc samples
t = linspace(0,T*numChirps*numCPI,N); % time axis, one frame
% t= 0:dt:dt*numADC*numChirps-dt;
t_onePulse = 0:dt:dt*numADC-dt;
numTX = 1;
numRX = 8;
Vmax = lambda/(T*4); % Max Unamb velocity m/s
DFmax = 1/2*PRF; % = Vmax/(c/fc/2); % Max Unamb Dopp Freq
dR = c/(2*BW); % range resol
Rmax = F*c/(2*slope); % TI's MIMO Radar doc
Rmax2 = c/2/PRF; % lecture 2.3
dV = lambda/(2*numChirps*T); % velocity resol, lambda/(2*framePeriod)
d_rx = lambda/2; % dist. between rxs
d_tx = 4*d_rx; % dist. between txs

N_Dopp = numChirps; % length of doppler FFT
N_range = numADC; % length of range FFT
N_azimuth = numTX*numRX;
R = 0:dR:Rmax-dR; % range axis
V = linspace(-Vmax, Vmax, numChirps); % Velocity axis
ang_ax = -90:90; % angle axis

%% Targets

r1_radial = 50;
tar1_angle = -15;
r1_y = cosd(tar1_angle)*r1_radial;
r1_x = sind(tar1_angle)*r1_radial;
v1_radial = 10; % velocity 1
v1_y = cosd(tar1_angle)*v1_radial;
v1_x = sind(tar1_angle)*v1_radial;
r1 = [r1_x r1_y 0];

r2_radial = 100;
tar2_angle = 10;
r2_y = cosd(tar2_angle)*r2_radial;
r2_x = sind(tar2_angle)*r2_radial;
v2_radial = -15; % velocity 2
v2_y = cosd(tar2_angle)*v2_radial;
v2_x = sind(tar2_angle)*v2_radial;
r2 = [r2_x r2_y 0];

tx_loc = cell(1,numTX);
for i = 1:numTX
   tx_loc{i} = [(i-1)*d_tx 0 0];
   scatter3(tx_loc{i}(1),tx_loc{i}(2),tx_loc{i}(3),'b','filled')
   hold on
end

rx_loc = cell(1,numRX);
for i = 1:numRX
   rx_loc{i} = [tx_loc{numTX}(1)+d_tx+(i-1)*d_rx 0 0];
   scatter3(rx_loc{i}(1),rx_loc{i}(2),rx_loc{i}(3),'r','filled')
end

tar1_loc = zeros(length(t),3);
tar2_loc = zeros(length(t),3);
tar1_loc(:,1) = r1(1) + v1_x*t;
tar2_loc(:,1) = r2(1) + v2_x*t;
tar1_loc(:,2) = r1(2) + v1_y*t;
tar2_loc(:,2) = r2(2) + v2_y*t;
% for i = 1:1000:length(t)
%     scatter3(tar1_loc(i,1),tar1_loc(i,2),tar1_loc(i,3),'m')
%     hold on
%     scatter3(tar2_loc(i,1),tar2_loc(i,2),tar2_loc(i,3),'k')
% end

%% TX

delays_tar1 = cell(numTX,numRX);
delays_tar2 = cell(numTX,numRX);
r1_at_t = cell(numTX,numRX);
r2_at_t = cell(numTX,numRX);
tar1_angles = cell(numTX,numRX);
tar2_angles = cell(numTX,numRX);
tar1_velocities = cell(numTX,numRX);
tar2_velocities = cell(numTX,numRX);
for i = 1:numTX
    for j = 1:numRX
        delays_tar1{i,j} = (vecnorm(tar1_loc-repmat(rx_loc{j},N,1),2,2)+vecnorm(tar1_loc-repmat(tx_loc{i},N,1),2,2))/c; 
        delays_tar2{i,j} = (vecnorm(tar2_loc-repmat(rx_loc{j},N,1),2,2)+vecnorm(tar2_loc-repmat(tx_loc{i},N,1),2,2))/c;
%         r1_at_t{i,j} = vecnorm(tar1_loc-repmat(rx_loc{j},N,1),2,2);
%         r2_at_t{i,j} = vecnorm(tar2_loc-repmat(rx_loc{j},N,1),2,2);
%         tar1_angles{i,j} = acosd(r1(2)./r1_at_t{i,j});
%         tar2_angles{i,j} = acosd(r2(2)./r2_at_t{i,j});
%         tar1_velocities{i,j} = v1 * cosd(tar1_angles{i,j});
%         tar2_velocities{i,j} = v2 * cosd(tar2_angles{i,j});
    end
end

%% Complex signal
phase = @(tx,fx) 2*pi*(fx.*tx+slope/2*tx.^2); % transmitted
phase2 = @(tx,fx,r,v) 2*pi*(2*fx*r/c+tx.*(2*fx*v/c + 2*slope*r/c)); % downconverted

% f_oneChirp =  slope*t(1:sum(t<=T));
% f_t = repmat(f_oneChirp,1,numChirps*numCPI)-(BW/2); % transmit freq
% f_t = BW/2*sawtooth(t/T*2*pi); 

fr1 = 2*r1(2)*slope/c; 
fr2 = 2*r2(2)*slope/c;
fd1 = 2*v1_radial*fc/c; % doppler freq
fd2 = 2*v2_radial*fc/c;
f_if1 = fr1 + fd1; % beat or IF freq
f_if2 = fr2 + fd2;


% mixed1 = cell(numTX,numRX);
% mixed2 = cell(numTX,numRX);
mixed = cell(numTX,numRX);
for i = 1:numTX
    for j = 1:numRX
        disp(['Processing Channel: ' num2str(j) '/' num2str(numRX)]);
        for k = 1:numChirps*numCPI
            phase_t = phase(t_onePulse,fc);
            phase_1 = phase(t_onePulse-delays_tar1{i,j}(k*numADC),fc); % received
            phase_2 = phase(t_onePulse-delays_tar2{i,j}(k*numADC),fc);
            
            signal_t((k-1)*numADC+1:k*numADC) = exp(1j*phase_t);
            signal_1((k-1)*numADC+1:k*numADC) = exp(1j*(phase_t - phase_1));
            signal_2((k-1)*numADC+1:k*numADC) = exp(1j*(phase_t - phase_2));
        end
        mixed{i,j} = signal_1 + signal_2;
    end
end

figure
subplot(3,1,1)
p1 = plot(t, real(signal_t));
title('TX')
xlim([0 0.1e-4])
xlabel('Time (sec)');
ylabel('Amplitude');
subplot(3,1,2)
p2 = plot(t, real(signal_1));
title('RX')
xlim([0 0.1e-4])
xlabel('Time (sec)');
ylabel('Amplitude');
subplot(3,1,3)
p3 = plot(t,real(mixed{i,j}));
title('Mixed')
xlim([0 0.1e-4])
xlabel('Time (sec)');
ylabel('Amplitude');

%% Post processing - 2-D FFT

RDC = reshape(cat(3,mixed{:}),numADC,numChirps*numCPI,numRX*numTX); % radar data cube
RDMs = zeros(numADC,numChirps,numTX*numRX,numCPI);
for i = 1:numCPI
    RD_frame = RDC(:,(i-1)*numChirps+1:i*numChirps,:);
    RDMs(:,:,:,i) = fftshift(fft2(RD_frame,N_range,N_Dopp),2);
end

figure
imagesc(V,R,20*log10(abs(RDMs(:,:,1,1))/max(max(abs(RDMs(:,:,1,1))))));
colormap(jet(256))
% set(gca,'YDir','normal')
clim = get(gca,'clim');
caxis([clim(1)/2 0])
xlabel('Velocity (m/s)');
ylabel('Range (m)');

%% CA-CFAR

numGuard = 2; % # of guard cells
numTrain = numGuard*2; % # of training cells
P_fa = 1e-5; % desired false alarm rate 
SNR_OFFSET = -5; % dB
RDM_dB = 10*log10(abs(RDMs(:,:,1,1))/max(max(abs(RDMs(:,:,1,1)))));

[RDM_mask, cfar_ranges, cfar_dopps, K] = ca_cfar(RDM_dB, numGuard, numTrain, P_fa, SNR_OFFSET);

figure
h=imagesc(V,R,RDM_mask);
xlabel('Velocity (m/s)')
ylabel('Range (m)')
title('CA-CFAR')

%% Angle Estimation - FFT

rangeFFT = fft(RDC(:,1:numChirps,:),N_range);

angleFFT = fftshift(fft(rangeFFT,length(ang_ax),3),3);
range_az = squeeze(sum(angleFFT,2)); % range-azimuth map

figure
colormap(jet)
imagesc(ang_ax,R,20*log10(abs(range_az)./max(abs(range_az(:))))); 
xlabel('Azimuth Angle')
ylabel('Range (m)')
title('FFT Range-Angle Map')
set(gca,'clim', [-35, 0])

doas = zeros(K,181); % direction of arrivals
figure
hold on; grid on;
for i = 1:K
    doas(i,:) = fftshift(fft(rangeFFT(cfar_ranges(i),cfar_dopps(i),:),181));
    plot(ang_ax,10*log10(abs(doas(i,:))))
end
xlabel('Azimuth Angle')
ylabel('dB')

%% Angle Estimation - MUSIC phased.Toolbox

% virt_Array = phased.ULA('NumElements',numRX*numTX,'ElementSpacing',d_rx);
% estimator = phased.MUSICEstimator('SensorArray',virt_Array,'OperatingFrequency',fc,...
%     'DOAOutputPort',true,'NumSignalsSource','Property','NumSignals',1,'ScanAngles',ang_ax);
% noise_pow = 3;
% M = numCPI; % # of snapshots
% 
% % sigma2 = 0.01; % Noise variance
% % n = sqrt(sigma2)*(randn(numRX*numTX,M) + 1j*randn(numRX*numTX,M))/sqrt(2);
% rangeFFT = fft(RDC,N_range);
% for i = 1:N_range
%     Rxx = zeros(numTX*numRX,numTX*numRX);
%     for m = 1:M
%        A = squeeze(sum(rangeFFT(i,(m-1)*numChirps+1:m*numChirps,:),2));
%        Rxx = Rxx + 1/M * (A*A');
%     end
% %     Rxx = Rxx + sqrt(noise_pow/2)*(randn(size(Rxx))+1j*randn(size(Rxx)));
%     [y, doas_music] = estimator(Rxx);
%     range_az_music(i,:) = y;
% end
% 
% figure
% colormap(jet)
% imagesc(ang_ax,R,20*log10(abs(range_az_music)./max(abs(range_az_music(:))))); 
% xlabel('Azimuth')
% ylabel('Range (m)')
% title('MUSIC Range-Angle Map')
% clim = get(gca,'clim');

%% Angle Estimation - MUSIC Pseudo Spectrum

d = 0.5;
M = numCPI; % # of snapshots

for k=1:length(ang_ax)
        a1(:,k)=exp(-1i*2*pi*(d*(0:numTX*numRX-1)'*sin(ang_ax(k).'*pi/180)));
end
    
for i = 1:K
    Rxx = zeros(numTX*numRX,numTX*numRX);
    for m = 1:M
       A = squeeze(RDMs(cfar_ranges(i),cfar_dopps(i),:,m));
       Rxx = Rxx + 1/M * (A*A');
    end

    [Q,D] = eig(Rxx); % Q: eigenvectors (columns), D: eigenvalues
    [D, I] = sort(diag(D),'descend');
    Q = Q(:,I); % Sort the eigenvectors to put signal eigenvectors first
    Qs = Q(:,1); % Get the signal eigenvectors
    Qn = Q(:,2:end); % Get the noise eigenvectors

    for k=1:length(ang_ax)
        music_spectrum(i,k)=(a1(:,k)'*a1(:,k))/(a1(:,k)'*(Qn*Qn')*a1(:,k));
    end
end

figure
hold on 
grid on
title('MUSIC Spectrum')
xlabel('Angle in degrees')
for k = 1:K
    plot(ang_ax,log10(abs(music_spectrum(k,:))));
end

%% Point Cloud

[~, I] = max(music_spectrum(2,:));
angle1 = ang_ax(I);
[~, I] = max(music_spectrum(1,:));
angle2 = ang_ax(I);

coor1 = [cfar_ranges(2)*sind(angle1) cfar_ranges(2)*cosd(angle1) 0];
coor2 = [cfar_ranges(1)*sind(angle2) cfar_ranges(1)*cosd(angle2) 0];
figure
hold on;
title('3D Coordinates of the targets')
scatter3(coor1(1),coor1(2),coor1(3),100,'m','filled','linewidth',9)
scatter3(coor2(1),coor2(2),coor2(3),100,'b','filled','linewidth',9)
xlabel('Range (m) X')
ylabel('Range (m) Y')
zlabel('Range (m) Z')
%% MUSIC Range-AoA map
rangeFFT = fft(RDC);
for i = 1:N_range
    Rxx = zeros(numTX*numRX,numTX*numRX);
    for m = 1:M
       A = squeeze(sum(rangeFFT(i,(m-1)*numChirps+1:m*numChirps,:),2));
       Rxx = Rxx + 1/M * (A*A');
    end
%     Rxx = Rxx + sqrt(noise_pow/2)*(randn(size(Rxx))+1j*randn(size(Rxx)));
    [Q,D] = eig(Rxx); % Q: eigenvectors (columns), D: eigenvalues
    [D, I] = sort(diag(D),'descend');
    Q = Q(:,I); % Sort the eigenvectors to put signal eigenvectors first
    Qs = Q(:,1); % Get the signal eigenvectors
    Qn = Q(:,2:end); % Get the noise eigenvectors

    for k=1:length(ang_ax)
        music_spectrum2(k)=(a1(:,k)'*a1(:,k))/(a1(:,k)'*(Qn*Qn')*a1(:,k));
    end
    
    range_az_music(i,:) = music_spectrum2;
end

figure
colormap(jet)
imagesc(ang_ax,R,20*log10(abs(range_az_music)./max(abs(range_az_music(:))))); 
xlabel('Azimuth')
ylabel('Range (m)')
title('MUSIC Range-Angle Map')
clim = get(gca,'clim');

%% Angle Estimation - Compressed Sensing

numTheta = length(ang_ax); % divide FOV into fine grid
B = a1; % steering vector matrix or dictionary, also called basis matrix
% s = ones(numTheta,1);
% psix = dftmtx(numTheta);
% inv_psix = conj(psix)/numTheta;
% cap_theta = B*psix; % random measurement of basis function
figure
hold on; grid on;
title('Angle Estimation with Compressed Sensing')
xlabel('Azimuth')
ylabel('dB')
for i = 1:K
    A = squeeze(RDMs(cfar_ranges(i),cfar_dopps(i),:,1));
    cvx_begin
        variable s(numTheta) complex; %alphax(numTheta,1) phix(numTX*numRX,numTheta)...
            %cap_theta(numTX*numRX,numTheta) %B(numTX*numRX,numTheta)%psix(numTheta,numTheta) %A(numRX*numTX,1) % A is the initial measurement
    %     cap_theta == phix * psix;
    %     minimize(norm(alphax,1))
    %     pow_p(norm(A-cap_theta*alphax,2),2) <= 1;
    %     norm(A-cap_theta*alphax,2) <= 1;

    %     minimize(norm(A-cap_theta*alphax,1))
        minimize(norm(s,1))
        norm(A-B*s,2) <= 1;
    cvx_end
    cvx_status
    cvx_optval

    plot(ang_ax,10*log10(abs(s)))
end

%% Compressed Sensing - Range AoA Map doesn't work

% for i = 1:N_range
%     A = squeeze(sum(rangeFFT(i,1:numChirps,:),2));
%     cvx_begin quiet
%         variable s(numTheta) complex; 
%         minimize(norm(s,1))
%         norm(A-B*s,2) <= 1;
%     cvx_end
%     range_az_cs(i,:) = s;
%     disp(['Processing Compressed Sensing ... ' num2str(i) '/' num2str(N_range)])
% end
% 
% figure
% colormap(jet)
% imagesc(ang_ax,R,20*log10(abs(range_az_cs)./max(abs(range_az_cs(:))))); 
% xlabel('Azimuth')
% ylabel('Range (m)')
% title('Compressed Sensing Range-Angle Map')
%% DOA estimation using Root SBA compress sensing
% figure; 
% hold on; grid on
% for i = 1:K
%     A = squeeze(RDMs(cfar_ranges(i),cfar_dopps(i),:,1));
%     [Pm_root,search_root]=Bayesian_DOA_root(A,ang_ax,M);
%     plot(search_root,Pm_root)
% end
