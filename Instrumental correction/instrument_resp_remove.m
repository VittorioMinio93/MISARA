%%%%%%%%%%%%% Instrument response correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nyquist: Nyquist frequency in Hz (1/2*sampling frequency)
% poles: instrument Poles from factory calibration worksheet (rad/s)
% zeroes: instrument Zeroes from factory calibration worksheet (rad/s) 
% bmv : coefficient for conversion from bits to Volts (V/counts) 
% vvel: coefficient for conversion from Volts to m/s (V s/m)
% k = calibration coefficient (rad/s)
% tracein: uncorrected seismic time series (seismic velocity)
% traceout: corrected seismic velocity trace

function traceout= instrument_resp_remove(nyquist,poles,zeroes,k,bmv,vvel,tracein)

% Before deconvolving instrument response we must get the seismogram into units of m/s
% There's digitiser conversion coefficient from bits to microVolts
% Then conversion from Volts to m/s

%Remove the best straight-fit line
 tracein = detrend(tracein);

%Deconvolving instrument response
if ~isnan(poles(1)) & ~isnan(zeroes(1)) 
[traceout] = pazrem(nyquist,poles,zeroes,k,tracein);
traceout=traceout';
else
traceout=tracein;
end
end

%==========================================================================
%%%%%%%%%%%%%%%%%%%%% Function to remove the instrument response %%%%%%%%%%
function [traceout] = pazrem(nyquist,poles,zeroes,k,tracein)
%Apply high pass filter and remove the best straight-fit line
n = 2;
Wn = [0.05/nyquist];
[d,c] = butter(n, Wn, 'high');
tracein = detrend(tracein);
tracein = filtfilt(d,c,tracein);
tracein = (detrend(tracein))';

%Gives the poles and zeros in Hz
[b,a] = zp2tf(zeroes,poles,k);

%Gives the (discrete) frequency response of the transfer function
j = nyquist/(length(tracein));
qq = [j(1):j:nyquist];
h =freqs(b,a,qq);

% PROBLEM WITH ABOVE : This section provides an unstable filter
% (This can be seen by the fact that three of the five poles are
% actually located outside of the z-plane's unit circle.
% Therefore the following is a routine to give a water-level type
% stabilisation.
% This requires that the amplitude spectrum does not fall below a
% given threshold BUT the phase spectrum can not change...
% The threshold chosen here (1.2e-2) is chosen due to the fact it is the
% abs(h) value at a period of 1000s

threshold = abs(h(end));

for xx = 1:length(h)

    if abs(h(xx)) < threshold
        imoverre = imag(h(xx))/real(h(xx));
        realh = sqrt((threshold*threshold)/(1+imoverre*imoverre));
        imagh = realh*imoverre;
        h(xx) = realh + i*imagh;
    end
end


% The phase spectrum of the above looks a bit suspect, however it is
% due to the fact that the different quadrants of the complex plane
% are not taken into account and therefore there is a change from
% -pi/2 to pi/2 as the function h crosses the imag axis. D.G 4/3/2003
% below take into the frequency domian apply the transfer function
% bring back to time domain.

matpazfreq = fft(tracein)./h;
matpaz = ifft(matpazfreq);
traceout = real(matpaz);

%Apply high pass filter and remove the best straight-fit line
n =2;
Wn = [0.05/nyquist];
[d,c] = butter(n, Wn, 'high');
traceout = detrend(traceout);
traceout = filtfilt(d,c,traceout);
traceout = detrend(traceout);
end