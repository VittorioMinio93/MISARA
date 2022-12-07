%%%%%%%%%%%%%Double cosine bell tapering%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x: signal vector
% npts: the analysis window for the MUSIC analysis (samples) 
% tb,te: beginning and ending taper in percent (integer)

function x=taper(x,npts,tb,te)
pio2=2*atan(1);
%Beginning tapering
if tb~=0
    n=round((npts*tb)/100);
    arg=pio2.*([0:n-1]./n);
    x(1:n)=x(1:n).*sin(arg);
end
%Ending tapering
if te~=0
    n=round((npts*te)/100);
    arg=pio2.*([0:n-1]./n);
    x(npts-n+1:npts)=x(npts-n+1:npts).*sin(flip(arg));
end
end