function FAP = pazToFreqResp(PAZ, omegas)
%PAZTOFREQRESP Convert Poles and Zeros to frequency response. 
% 
% FAP = PAZTOFREQRESP(PAZ,omegas) Convert poles and
% zeros (in units of rad/s) to frequency response (frequency/amplitude pairs) at angular freqs
% specified by vector omegas (rad/s).

% Author: Silvio De Angelis, School of Environmental Sciences, University
% of Liverpool
% $Date: 2011-09-17 00:00:00 $
% $Revision: 1.0 $

%INITIALIZE FAP STRUCTURE AND POPULATE FIELDS
FAP.frequencies = reshape(omegas/(2*pi),1,numel(omegas));
FAP.values = [];
FAP.freqUnits = 'Hz';
% NORMALIZE POLES AND ZEROS
%PAZ.normalization = 1./abs(polyval(poly(PAZ.zeros),2*pi*1i)/polyval(poly(PAZ.poles),2*pi*1i));

% CALCULATE FREQUENCY RESPONSE AT COMPLEX omegas
[a,b] = zp2tf(PAZ.zeros,PAZ.poles,PAZ.normalization)
FAP.values = freqs(a,b,omegas);
FAP.phase = (angle(FAP.values))*180/pi;





