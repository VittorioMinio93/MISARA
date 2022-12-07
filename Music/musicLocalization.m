%%%%%%%%%%%%%  MUSIC algorithm  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dat: signals matrix
% ndat:length of the signals
% x: station coordinates matrix
% n, nsta: number of station
% dt, delta_time: sampling interval
% fsam: sampling rate 
% freq0: the central frequency of the band of interest (Hz)
% nfreq: the number of frequencies of the band of interest
% fstep: frequency analyisis step
% f0: focusing frequency
% ifo: index of the focusing frequency
% cf: index of frequency analysis
% m: the hamming value
% R: the R value
% slowx0,slowy0,pmax: the limits of the grid slowness (s/km) 
% ds, pinc: the size step of the grid slowness (s/km) 
% ns: dimensions of slowness grid
% ninicio: the time origin of the analysis window (s)  
% tprime,lwin: the analysis window for the MUSIC analysis (s) 
% num_time,nwin: number of analysis windows
% ntime: number of time window reference
% advance: the advancement of the analysis window for the MUSIC analysis (%) 
% time_step: temporal step of advance
% Nsig: the number of sources computed 
% eigperc: the threshold value of the sources selection
% freq_rms: frequency interval analysis for rms analysis
% ta, tb: begining and ending taper in percent (integer)
% a,aa: haming value
% cu: spectrums matrix
% df: frequency step of spectrums
% tt: the focusing matrix
% amp, phi: real and imaginary parts of the spectrums
% rr: the focused spatial matrix for a specific frequency
% D,V: eigenvalues and eigenvectors matrixs
% L: projection length
% E: eigenvectors selected matrix
% dataflag: active/deactive normalization of the signals to PGA
% flagla:  Active/Deactive the normalization of the spectrums 
% flaglb: Active/Deactive the normalization of the focused cross-spectral matrix 
% flagsig: method to determine ranges of eigenvectors to use
% Nsig1,NSig2: number of signal sources/ signal noise
% BAZI,bazi: backazimuth vector/matrix  (degrees)
% SLOW,slow: ray parameter vector/matrix (s/km) 
% power: fk-power spectrum 
% peak: main peaks of the fk-power spectrum
% FKPOW,fkpower: fk power vector/matrix 
% INC,Inc: incidence vector/matrix (degrees)
% SXMAX/sxmax, SYMAX/symax: horizontal slowness vectors/matrixs (s/km) 
% Time, tiempo: time vectors
% TMP, TMP_perc: the cumulative matrix/vector 
% iday,ihr,imin,ssec: time references 

function [Time,SLOW,BAZI,FKPOWER,TMP,SXMAX,SYMAX,freq_rms]=musicLocalization(dat,x,ninicio,fsam,dt,lwin,advance,R,nfreq,freq0,m,pmax,pinc,Nsig1,Nsig2,nsta,a,tt,eigperc)
row1=30;
row=2*30;
%Set temporal variables
iday=0;
ihr=0;
imin=0;
ssec=ninicio*dt;
%Set the analysis parameters
time_step=lwin*advance/fsam;ndat=length(dat);
nwin=round((ndat)/(time_step*fsam));
num_time=nwin;
tprime=lwin;
tb=10;te=10;
flag1a='N';flag1b='N';
fstep=fsam/tprime;
ifo=fix(nfreq/2)+1;
flagsig='N';
slowx0=-pmax;
slowy0=-pmax;
ds=pinc;
ns=2*round(pmax/pinc)+1;
%Initialize initial slowness
sx_steer=0;algn_sx=0;algn_sy=0;
sy_steer=0;
sx_shift=algn_sx+sx_steer;
sy_shift=algn_sy+sy_steer; 
tmin=0;
dataflag=0;
n=nsta;
Time=NaN(nwin,1);%Initialize the time vector
SLOW=NaN(nwin,nsta);%Initialize the slowness vector
BAZI=NaN(nwin,nsta);%Initialize the backazimuth vector
FKPOWER=NaN(nwin,nsta);%Initialize the fk-power vector
TMP=NaN(nwin,nsta);%Initialize the cumulative vector
SXMAX=NaN(nwin,nsta);SYMAX=NaN(nwin,nsta);%Initialize the horizontal slowness vector
%Set the focusing frequency
df=1./(tprime*dt);
[fo,~]=setfreq(ifo,freq0,df,fstep);
% Form spatial covariance matrix of DFT
if (nfreq==1) | (ifo==0)
    aa=1./[1:nfreq];
else
    nn=nfreq-ifo;
    if ifo>nn;nn=ifo-1;end 
    [aa,~]=haming(nn);   
    aa=aa./sum(aa);
end
% Move the analysis window through the seismic traces
ntime=0; %for the next window
nskip=ninicio+1; %for the next window
while ntime<num_time & nskip+tprime-1<=ndat
ntime=ntime+1; %for the next window
%%Initialize focused covariance matrix
s=zeros(2*n,2*n);
%%Read spectrums
[df,cu,~]=rdf(dat,nsta,dt,nskip,tprime,tmin,tb,te,dataflag);
delta_time=dt;freq_rms=zeros(nfreq,1);
% Loop through each frequency bin
for j=1:nfreq
    %Set frequency
    [frq,cf]=setfreq(j,freq0,df,fstep);
      %Return real and imaginary parts of the spectrums
      [phi,amp]=savspt(cf,m,flag1a,n,cu);
      %Estimate of the focused spatial matrix for current frequency
      rr=crspt (n,m,phi,amp,a,tt);
      %Update the weighted sum of the focused spatial matrix 
      jj = nn + 1 + j - ifo;
      s=s+rr.*aa(jj);
      freq_rms(j)=frq;
end
%Normalize focused cross-spectral matrix 
if (flag1b=='Y') | (flag1b=='y') 
    for i=1:n
        denom=s(i,i);
        for j=1:n
            s(i,j)=s(i,j)/denom;
            s(i,j+n)=s(i,j+n)/denom;
            s(i+n,j)=s(i+n,j)/denom;
            s(i+n,j+n)=s(i+n,j+n)/denom;
        end
    end
end
%Prewhittening
for i=1:2*n
    s(i,i) = s(i,i) + R;    
end
%Eigen-decomposition 
[V,Dd]=eig(s);
D=diag(Dd)';
%Separate eigenvectors into signals and noise subspace
[E,L,Nsig,TMP_perc,Nsig1,Nsig2]=signoise(row,row1,D,V,fo,flagsig,n,nfreq,tprime,Nsig1,Nsig2,eigperc);
%Evalutate reciprocal of the projection length over a grid slowness
power=eq12 (n,fo,E,L,x,row1,slowx0,slowy0,ds,ns,flagsig,sx_shift,sy_shift,Nsig,Nsig1,Nsig2);
%Normalize fk-power spectrum
power=power./max(power(:));
%Find the main peaks in the fk-power spectrum
peak=fpeak(ns,ds,slowx0,slowy0,Nsig,power);
%Return the Output of the analysis
[tiempo,slow,bazi,fkpower,sxmax,symax]=output(delta_time,lwin,imin,ssec,Nsig,peak,ntime,fo,x,n,E);
%Concatenate the results
Time(ntime)=tiempo;
SLOW(ntime,1:Nsig)=slow;
BAZI(ntime,1:Nsig)=bazi;
FKPOWER(ntime,1:Nsig)=fkpower;
TMP(ntime,:)=TMP_perc;
SXMAX(ntime,1:Nsig)=sxmax;
SYMAX(ntime,1:Nsig)=symax;
%%Increment for the next window analysis
    nskip=nskip+round(time_step/delta_time);
    ssec=ssec+time_step;
    if ssec>=60
       ssec=ssec-60;
       imin=imin+1;
       if imin>=60
           imin=imin-60;
           ihr=ihr+1;
           if ihr>=24
               ihr=ihr-24;
                iday=iday+1;
           end 
       end  
    end  
end
end