function trace = detrend(trace)
% DETREND  detrend seismogram
%
%   TRACE= DETREND(TRACE) detrend signal by subtracting a line through the first
%     and last point of the trace

% Author: Silvio De Angelis, School of Environmental Sciences, University
% of Liverpool
% $Date: 2011-09-17 00:00:00 $
% $Revision: 1.0 $
trace = trace(:)';
ndat = length(trace);
extr = [trace(1) trace(end)];
line = (extr(1)+(1:ndat)*(extr(2) - extr(1)) / (ndat - 1) );
line = line(:)';
trace = trace-line;