%%%%%%%%%%%%% Set the frequency analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% freq0: the central frequency of the band of interest (Hz)
% fstep: frequency analyisis step
% df: frequency step of spectrums
% j: frequency index respect on freq0
% freq: frequency analysis
% cf: index of frequency analysis

function [freq,cf]=setfreq(j,freq0,df,fstep)
freq = freq0 + (j-1)*fstep;
cf =round(freq/df+0.5);
freq = (cf-1)*df;
end