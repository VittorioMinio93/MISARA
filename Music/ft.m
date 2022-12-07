%%%%%%%%%%%%%Fast Fourier Transform algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u: signal vector
% dt: sampling interval
% npts: the analysis window for the MUSIC analysis (samples) 
% tb,te: begining and ending taper in percent (integer)
% dataflag: active/deactive normalization of the signals to PGA
% cu1: complex signal
% c: spectrums matrix
% df: frequency step of spectrums

function [c,df]=ft(u,npts,nmin,dt,tb,te,dataflag)
         %Rescale
         u=u-mean(u);
         %Tapering
         u=taper(u,npts,tb,te);
         %Active/Deactive normalization of the signals to PGA
         accmax=0;
         if dataflag<0
             for i=1:npts
                 if (abs(u(i))>accmax);accmax(i)=abs(u(i));end
             end
             u=u./accmax; 
         end 
         % Apply FFT algorithm
         cu1=complex(u);
         cu1=fft(cu1,npts);
         df=1./(npts*dt);
         c=cu1.*dt;      
end