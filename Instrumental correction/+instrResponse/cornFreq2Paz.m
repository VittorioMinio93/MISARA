function PAZ = cornFreq2Paz(fc, damp)
%CORNFREQ2PAZ  Convert corner frequency and damping to poles and zeros (in units of rad/s). 2 zeros at
%   position (0j, 0j) are given as output (m/s).

% Author: Silvio De Angelis, School of Environmental Sciences, University
% of Liverpool
% $Date: 2011-09-17 00:00:00 $
% $Revision: 1.0 $

PAZ.poles = [-(damp+sqrt(1-damp^2)*1j)*2*pi*fc];
PAZ.zeros = [0j,0j];
PAZ.normalization = 1/abs(polyval(poly(PAZ.zeros),2*pi*1i)/polyval(poly(PAZ.poles),2*pi*1i));