%%%%%%% Calculation of the delay time through the Cross-Correlation%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y1,y2: pair of signals
% splineflag: Activation/Deactivation cubic spline interpolation value (0 o 1)
% dt: sampling interval
% maxlags: maxlags: the max lags for the Correlation coefficient calculation (s)
% lags: temporal lags
% xc: Cross-Correlation coefficients 
% x: lags interval
% y: Cross-Correlation coefficients interval
% xx: resampled lags  interval
% yy: interpolated Cross-Correlation coefficients interval
% xcmax: maximum cross-Correlation coefficient
% delay: delay time values at the maximum cross-Correlation coefficient

function [delay xcmax]=GetDelay(y1,y2,splineflag,dt,maxlags)
maxlags=round(maxlags/dt);%Expressed in samples 
splineinterval=10;
splinefactor=10;
x=zeros(2*splineinterval+1,1);%Initialize the vector of the lags interval
y=zeros(2*splineinterval+1,1);%Initialize the vector of the  cross-correlation coefficients interval 
nx=length(x);
xx=zeros(nx+(nx-1)*splinefactor,1);%Initialize the vector of the resampled lags
yy=zeros(nx+(nx-1)*splinefactor,1);%Initialize the vector of the interpolated Cross-Correlation coefficients
%Calculation of the Cross-Correlation coefficients
[xc lags]=xcorr(y1,y2,maxlags,'coef');
[xcmax imax]=max(xc);%Estimation of the maximum value

%Activation/Deactivation of the cubic spline interpolation 
if splineflag==1 
    i1=imax-splineinterval; i2=imax+splineinterval;%Limits of the spline interval 
    %Control limits
    if i1 < 1; i1=1; end
    if i2 > length(lags); i2=length(lags); end
    %Apply the cubic spline interpolation
    x=(lags(i1:i2))*dt; 
    xx=linspace(min(x),max(x),splinefactor*nx+1);
    y=xc(i1:i2);
    yy=spline(x,y,xx);
    %Calculation of the delay time value at the maximum cross-Correlation coefficient interpolated 
    [tmax imax]=max(yy);
    delay=xx(imax);
else
    %Calculation of the delay time value at the maximum cross-Correlation coefficient 
    delay=lags(imax)*dt;
end
end