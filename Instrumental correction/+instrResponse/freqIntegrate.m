function traceOut = freqIntegrate(traceIn,fs)
%FREQINTEGRATE perform integration in the frequency domain 

%INITIALIZE FFT
traceOutFFT = zeros(1,length(traceIn));

yfft = fft(traceIn);

for j=1:length(traceIn)
    omega = 2*pi*(fs/2/length(traceIn))*j;
    traceOutFFT(j) = imag(yfft(j)) / omega + real(yfft(j)) / (omega * 1i);
end

traceOut = ifft(traceOutFFT,'symmetric');
