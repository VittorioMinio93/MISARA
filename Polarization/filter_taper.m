%%%%%%%%%%%%%Application of the polarization analysis%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% az: azimuth vectors (degrees)
% incid: incidence vectors (degrees)
% rect: rectilinearity vectors
% sis: three components matrix 
% wl, ws: windows analysis (s)
% f1, f2: frequency analysis (Hz) 
% fs: sampling rate
% nyq: nyquist number
% nsamp: number of samples

function [az,incid,rect,LL,TL] = filter_taper (sis,fs,f1,f2,wl,ws)
%Rescale the seismic traces
%Loop through the three components 
for k=1:3
    sis(:,k)=sis(:,k)-mean(sis(:,k));
end
%Tapering
[nsamp col]=size(sis);%dimensions of the 3D trace matrix
kk=round(nsamp*1/50); % 2% of the trace 
%Loop through the edges (2%) of the samples
for i=1:kk
    sis(i,:)=sis(i,:)*(1-cos(pi*(i-1)/(2*(kk-1))));
    sis(nsamp+1-i,:)=sis(nsamp+1-i,:)*(1-cos(pi*(i-1)/(2*(kk-1))));
end

%Filter the seismic traces with a band-pass
nyq=fs/2;
w1=f1/nyq; w2=f2/nyq;
W=[w1 w2];
[a,b]=butter(2,W);
sis=filtfilt(a,b,sis);
%Recall the Jurkevics algorithm function
[az,incid,rect,LL,TL] = polarization(sis,fs,wl,ws);
