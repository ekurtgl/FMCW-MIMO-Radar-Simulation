function [RDM_mask, cfar_ranges, cfar_dopps, K] = ca_cfar(RDM_dB, numGuard, numTrain, P_fa, SNR_OFFSET)
        % e.g. numGuard =2, numTrain =2*numGuard, P_fa =1e-5, SNR_OFFSET = -15
numTrain2D = numTrain*numTrain - numGuard*numGuard;
RDM_mask = zeros(size(RDM_dB));

for r = numTrain + numGuard + 1 : size(RDM_mask,1) - (numTrain + numGuard)
    for d = numTrain + numGuard + 1 : size(RDM_mask,2) - (numTrain + numGuard)
        
        Pn = ( sum(sum(RDM_dB(r-(numTrain+numGuard):r+(numTrain+numGuard),d-(numTrain+numGuard):d+(numTrain+numGuard)))) - ...
            sum(sum(RDM_dB(r-numGuard:r+numGuard,d-numGuard:d+numGuard))) ) / numTrain2D; % noise level
        a = numTrain2D*(P_fa^(-1/numTrain2D)-1); % scaling factor of T = α*Pn
        threshold = a*Pn;
        if (RDM_dB(r,d) > threshold) && (RDM_dB(r,d) > SNR_OFFSET)
            RDM_mask(r,d) = 1;
        end
    end
end

% figure(2)
% imagesc(RDM_mask);
% title('CA-CFAR')

[cfar_ranges, cfar_dopps]= find(RDM_mask); % cfar detected range bins

%% remaining part is for target location estimation
rem_range = zeros(1,length(cfar_ranges));
rem_dopp = zeros(1,length(cfar_dopps));
for i = 2:length(cfar_ranges)
   if (abs(cfar_ranges(i) - cfar_ranges(i-1)) <= 5) && (abs(cfar_dopps(i) - cfar_dopps(i-1)) <= 5)
       rem_range(i) = i; % redundant range indices to be removed
       rem_dopp(i) = i; % redundant doppler indices to be removed
   end
end
rem_range = nonzeros(rem_range); % filter zeros
rem_dopp = nonzeros(rem_dopp); % filter zeros
cfar_ranges(rem_range) = [];
cfar_dopps(rem_dopp) = [];
K = length(cfar_dopps); % # of detected targets
end