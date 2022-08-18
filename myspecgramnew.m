function [out1] = myspecgramnew(data,window,nfft,shift)

N = floor((length(data)-window-1)/shift); % 104
for i=1:N
    
    tmp =(fft(data((i-1)*shift+1:(i-1)*shift+window).'.*hann(window),nfft));
    out1(:,i) = tmp;
end
