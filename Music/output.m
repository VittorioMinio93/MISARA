%%%%%%%%%%%%%%Return the results of the analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dtime: sampling interval
% lwin: the analysis window for the MUSIC analysis (s) 
% imin: time minutes of  reference
% ssec: time seconds of  reference
% Nsig: the number of sources computed 
% peak: main peaks of the fk-power spectrum
% ntime: number of time window reference 
% freq: frequency analysis
% x: station coordinates matrix
% n: numeber of stations
% e: eigenvectors selected matrix
% tiempo: time vector
% slow: ray parameter matrix (s/km) 
% baz,bazi: backazimuth matrix  (degrees)
% fkpower: fk power matrix
% sxmax/Sxmax,symax/Symax: horizontal slowness vectors/matrixs (s/km) 
% tphase: theorical phase
% ophase: observed phase
% diff: difference betweem the i-th theorical phase and the observed phase
% re: real part
% im: imaginary part 

function [tiempo,slow,bazi,fkpower,Sxmax,Symax]=output (dtime,lwin,imin,ssec,Nsig,peak,ntime,freq,x,n,e)
tiempo=imin*60+ssec+lwin*dtime;
for i=1:Nsig
    fkpower(i)=peak(1,i);
    sxmax=peak(2,i);
    symax=peak(3,i);
    slow(i)=sqrt(sxmax*sxmax+symax*symax);
    baz=90-atan2(symax,sxmax).*180./pi;
    k=find(baz<0); baz(k)=baz(k)+360;
    bazi(i)=baz;
    Sxmax(i)=sxmax;
    Symax(i)=symax;
end
%Determine miss-fit if Nsig = 1
if Nsig==1
%     pi=atan(1.0)*4.0;
    ome=2.0*pi*freq;
    j=1;
    %calculate theoretical phase delay 'tphase' relative to the reference
    for i=1:n
        tphase(i)=ome*(peak(2,j)*x(1,i)+peak(3,j)*x(2,i));
    end
    refph = tphase(1);
    for i=1:n
        tphase(i)=tphase(i)-refph;
    end
    %calculate observed phase delay and use the theoretical phase delay to unwrapp the phase
    for i=1:n
      re = real(e(i,j));
      im = imag(e(i,j));

    if re==0
        if im<0;ophase(i)=pi * 3.0 /2.0;end
        if im>0;ophase(i)=pi/2.0;end
    else
     if (im>0) & (re>0);ophase(i)=atan2(im,re);end 
     if (im>0) & (re<0);ophase(i)=pi-atan2(im,-1.0*re);end 
     if (im<0) & (re<0);ophase(i)=atan2(-1.0*im,-1.0*re)-pi;end 
     if (im<0) & (re>0);ophase(i)=-1.0*atan2(-1.0*im,re);end     
    end
    end
    refph = ophase(1);
    diff=2*pi;ip=0;
    for i=1:n
        if abs(diff)>=2*pi
        ipi=0;
        ophase(i)=refph-ophase(i);
        end
        tmp = ophase(i)+ipi*2.0*pi*sign(tphase(i));
        diff= tmp-tphase(i);
        if abs(diff)<2*pi
            ophase(i)=tmp;
        else
            ip=ip+1;
        end
    end    
end
end