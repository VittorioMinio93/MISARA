%%%%%%%%%%%%%Particle motion algorithm%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path: the path of the coordinate file (.mat) 
% staref: the name of the reference station 
% staidx: station index of reference
% station: station names
% xs,ys,zs,utm_zone:  station coordinates in UTM system
% LAT,LON,ELE: station coordinates in wgs84 system
% cooord: coordinates of the station of reference
% sis: the three components matrix
% dt: sampling interval
% fs: sampling rate
% nyq: nyquist number
% f1, f2: frequencies of the band-pass filter
% wls: analysis window (s) 
% coord_sys: the coordinate system (metrics/degrees)
% text1: Window message
% X_tot, Y_tot, Z_tot: Particle motions matrixs

function [X_tot,Y_tot,Z_tot,cooord]=ParticleMotion(path,sis,fs,f1,f2,wls,staref,coord_sys,text1)
%Initialize the station coordinates vectors
LAT=[];LON=[];ELE=[];xs=[];ys=[];zs=[];utm_zone=[];
%Initialize the coordinates vector of the reference station 
cooord=[];
%Load the stations coordinates
load(path)
%Set the station index of reference 
staidx=strcmp(station,staref);
%Set the coordinates vector of the reference station on the basis of the geographic coordinates system 
switch coord_sys
    case 'm'
        %Error control
        if isempty(xs);set(text1, 'String','Error. No correct coordinate file.');drawnow;return;end
         coords=[xs';ys';zs'];
         cooord=[coords(1, staidx),coords(2, staidx),coords(3, staidx)];
    case 'degrees'
        %Error control
        if isempty(LAT);set(text1, 'String','Error. No correct coordinate file.');drawnow;return;end
        cooord=[LAT(staidx),LON(staidx),ELE(staidx)];
        %Convert WGS84 coordinates (Latitude, Longitude) into UTM coordinates
        [xs,ys,~] = deg2utm(LAT,LON);
        %Express coordinates in km 
        xs=xs./1000;
        ys=ys./1000;
        zs=ELE./1000;
        coords=[xs';ys';zs'];        
end
%Rescale the seismic traces
%Loop through the three components 
for k=1:3
    sis(:,k)=sis(:,k)-mean(sis(:,k));
end
%Tapering
[nsamp col]=size(sis);%dimensions of the 3D trace matrix
kk=round(nsamp*1/50);% 2% of the trace 
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
sis=filter(a,b,sis);

%Set the analysis parameters
N=size(sis,1);wl=round(wls*fs);
i1=0; i2=0;k=0;
X_tot=[];Y_tot=[];Z_tot=[];%Initialize the Particle motion matrixs
%Move the analysis window through the seismic traces
while i2 < N-wl || i2== N-wl
    k=k+1;
    i1=(1+(k-1)*wl); i2=i1+wl-1;%Limits of the analysis window
    sisa=sis(i1:i2,:);%Cutted the seismic traces
    np=length(sisa);
    %Calculation of the Partcle motions
    sisam=sisa./repmat(max(sisa),np,1);
    X=sisam(:,3)+coords(1, staidx);
    Y=sisam(:,2)+coords(2, staidx);
    Z=sisam(:,1)+coords(3, staidx);
    %Eventually, convert UTM coordinates into WGS84 coordinates (Latitude, Longitude) 
    if coord_sys=='degrees'
        [Y,X] = utm2deg(X.*1000,Y.*1000,repmat(utm_zone(1,:),size(Z,1),1));
        %Expressed in meters
        Z=Z.*1000;
    end
    %Concatenate the results
    X_tot=[X_tot X];
    Y_tot=[Y_tot Y];
    Z_tot=[Z_tot Z];
end
end