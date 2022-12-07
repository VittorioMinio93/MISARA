%%%%%%%%%%%%%Calculation of the spectrums %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dat: signals matrix
% u: signal cutted
% N: number of station
% dt: sampling interval
% nskip: temporal advancing
% T,Tprime: the analysis window for the MUSIC analysis (samples) 
% tb,te: begining and ending taper in percent (integer)
% dataflag: active/deactive normalization of the signals to PGA
% c,cu: spectrums matrix
% df: frequency step of spectrums

function [df,cu,tdelay]=rdf(dat,N,dt,nskip,Tprime,Tmin,tb,te,dataflag)
         nfirst=1;
         for i=1:N
             ista=i;
             tdelay(ista)=0;%Initial temporal delay
             u=dat(ista,nskip:nskip+Tprime-1);
             T=Tprime;
             %Apply the fast fourier transform
             [c,df]=ft(u,T,Tmin,dt,tb,te,dataflag);
             cu(:,ista)=c';
         end        
end