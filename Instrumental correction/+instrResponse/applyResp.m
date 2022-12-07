function traceOut = applyResp(traceIn, nyq, filterBand, PAZ)
%
%APPLYRESP Apply instrument response to seismogram

% Author: Silvio De Angelis, School of Environmental Sciences, University
% of Liverpool
% $Date: 2011-08-17 00:00:00 $
% $Revision: 1.0 $

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

%IF NO PAZ IS GIVEN PAZ == WOOD-ANDERSON DISPLACEMENT RESPONSE.
if nargin < 4
    
    %WOOD-ANDERSON RESPONSE PAZ (NOTE THIS IS A DISPLACEMENT RESPONSE,
    % INPUT TRACE SHOULD BE DISPLACEMENT)
    PAZ.poles = [-5.49779-5.60886*1i; -5.49779+5.60886*1i];
    PAZ.zeros = [0+0*1i;0+0*1i];
    PAZ.gain  = 2080;
    PAZ.sensitivity = 1;
    PAZ.normalization = 1/abs(polyval(poly(PAZ.zeros),2*pi*1i)/polyval(poly(PAZ.poles),2*pi*1i));
    
end

%DETREND DATA
traceIn = traceIn(:)';
traceIn = instrResponse.detrend(traceIn);
dataLength = numel(traceIn);


%CREATE TAPER FILTER
taperWin  = instrResponse.cosTaper(dataLength,0.1);

%TAPER AND ZERO-PAD DATA TO INCREASE THE NUMBER OF SAMPLES X3
traceIn = traceIn.*taperWin(:)';
traceIn = [ zeros(1,dataLength) traceIn zeros(1,dataLength) ];
dataLength = numel(traceIn);

%DEAL WITH EVEN AND ODD NUMBER OF SAMPLES SEPARATELY
PAZ.poles = PAZ.poles(:);
PAZ.zeros = PAZ.zeros(:);


if floor(dataLength/2) == ceil(dataLength/2)
    
    %GET [O NYQ] FREQUENCIES IN RAD/SEC
    halfDataLength = dataLength/2;
    omega = 2*pi*[0:(halfDataLength)]*(1/(halfDataLength))*nyq;
    
    %GET FREQUENCY-AMPLITUDE PAIRS (FAP) FROM PAZ
    FAP = instrResponse.pazToFreqResp(PAZ,omega);
    
    %GET RESPONSE FOR NEGATIVE FREQUENCIES. THIS IS NEEDED LATER AS DECONVOLUTION IS
    %PERFORMED IN FREQUENCY DOMAIN WHERE FFT(DATA) IS TWO-SIDED.
    omega2 = [ -omega(halfDataLength+1:-1:2) omega(1:halfDataLength) ];
    
    %GET INVERSE RESPONSE FOR BOTH POSITIVE AND NEGATIVE FREQUENCIES AROUND
    %ZERO-OFFSET FREQUENCY.
    respFull = ([ real(FAP.values(end:-1:2)) real(FAP.values(1:end-1)) ] + 1i*[ -imag(FAP.values(end:-1:2)) imag(FAP.values(1:end-1)) ]);
    
else
    
    %THIS BIT IS THE SAME AS BEFORE BUT FOR ODD DATA LENGTH
    halfDataLength = (dataLength+1)/2;
    omega = 2*pi * [0:(halfDataLength-1)] * (2/dataLength) * nyq;
    FAP = instrResponse.pazToFreqResp(PAZ,omega);
    
    omega2 = [ -omega(end:-1:2) omega(1:end) ];
    respFull = ([ real(FAP.values(end:-1:2)) real(FAP.values(1:end)) ] + 1i*[ -imag(FAP.values(end:-1:2)) imag(FAP.values(1:end)) ]);
    
end

%DECONVOLVE INSTRUMENT RESPONSE BY DIVISION IN THE FREQUENC DOMAIN (MULTIPLICATION
% OF DATA BY INVERSE RESPONSE)
traceOut = PAZ.gain*(real(ifft(ifftshift(fftshift(fft(traceIn)).*(respFull)))));

%REMOVE INITIAL ZERO-PADDING
traceOut = traceOut( (dataLength/3)+1:2*(dataLength/3) );

%BANDPASS FILTER DATA TO REMOVE VERY LONG-PERIOD SIGNAL
[d,c] = butter(1, filterBand/nyq, 'bandpass');
traceOut = filtfilt(d,c,traceOut);
