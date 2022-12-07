function tapWin = cosTaper(npts, p)
%COSTAPER Cosine taper.
%   TAP = COSTAPER(NPTS, P) returns an NPTS-point cosine tapered window in a column vector, TAP.
%   The p parameter specifies the ratio of taper to constant sections. This ratio
%   is normalized to 1 (i.e., 0 < R < 1). If omitted, p is set to 0.10.

% Author: Silvio De Angelis, School of Environmental Sciences, University
% of Liverpool
% $Date: 2011-09-17 00:00:00 $
% $Revision: 1.0 $

error(nargchk(1,2,nargin,'struct'));

% Default value for R parameter.
if nargin < 2 || isempty(p),
    p = 0.100;
end

if  p <= 0,
    tapWin = ones(n,1);
elseif p >= 1,
    tapWin = hann(n);
else
    t = linspace(0,1,npts)';    
    per = p/2;
    tl = floor(per*(npts-1))+1;
    th = npts-tl+1;
    % Window is defined in three sections: taper, constant, taper
    tapWin = [ ((1+cos(pi/per*(t(1:tl) - per)))/2);  ones(th-tl-1,1); ((1+cos(pi/per*(t(th:end) - 1 + per)))/2)];
end

