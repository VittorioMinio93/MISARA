%%%%%%%%%%%%% JackKnife method for discrete analysis%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Z: signals matrix
% splint: the Spline interpolation value (1 o 0) 
% dt: sampling interval
% maxlags: the max lags for the Correlation coefficient calculation (samples)
% X: station coordinates matrix
% % wl,ws: windows analysis (samples)
% aaz: azimuth vector
% baz,bbaz: backazimuth vectors 
% rayp,rrap: ray paramater vectors
% sx,ssx: horizontal slowness vectors
% sy,ssy: horizontal slowness vectors
% NStation: number of station used
% Z_1: signals matrix leaving out one station at time 
% X_1: station coordinates matrix leaving out one station at time 
% delay: delay time 
% xcmax,Xcmax: cross-correlation coefficient vectors 
% params, par,PARR: parameters computed by omitting one station + parameters computed with all stations
% J: pseudovalues
% J_valor: JackKnife estimator of the parameters computed by omitting one station 
% J_diff: squared difference between the pseudovalues and the JackKnife estimator 
% J_error: standard error of the JackKnife estimator

function J_error=JackKnife_low(Z,splint,dt,maxlags,X,wl,ws,baz,rayp,sx,sy)
NStation=size(Z,2);
%%Loop through the number of stations
for gg=NStation:-1:1
    %Leave out the gg-th station 
    Z_1=Z;
    Z_1(:,gg)=[];
    X_2=X;
    X_2(gg,:)=[];
    ip=0;N=length(Z_1);nfile=size(Z_1,2);%Define the new number of stations (NStation-1)
    %Loop through all independent pairs of stations-1 
    for is1=1:nfile-1
        for is2=is1+1:nfile
            ip=ip+1;
            i1=0; i2=0;k=0;  
            %Move the analysis window through the seismic traces
            while i2 < N-wl || i2== N-wl
                k=k+1;
                i1=(1+(k-1)*ws); i2=i1+wl-1;%Limits of the analysis window
                y1=Z_1(i1:i2,is1).*tukeywin(wl); y2=Z_1(i1:i2,is2).*tukeywin(wl);%tapered cosine windows
                %Calculation of the Cross-correlation coefficients and the delay times 
                [delay(ip,k) xcmax(ip,k)]=GetDelay(y1,y2,splint,dt,maxlags);
            end
        end
    end
    %Calculation of the kinematic parameters and time vector
    [ssx ssy aaz rrayp]=FitPlane(X_2,-delay);
    Xcmax=mean(xcmax)';% Meaning all pairs of coefficients
    %Selection of discrete values
    [Xval,imax]=max(Xcmax);
    ssx=ssx(imax);ssy=ssy(imax);aaz=aaz(imax);rrayp=rrayp(imax);
    %Convert azimuth in backazimuth
    bbaz=aaz+180;
    %Concatenate the new parameters computed by omitting one station at time
    params(gg,:)=[{bbaz},{rrayp},{ssx},{ssy}];
end
%Concatenate the parameters computed with all stations
params(NStation+1,:)=[{baz},{rayp},{sx},{sy}];
npar=length(baz);%length of ZLC vectors
%Loop through the samples of the parameters
for jj=1:npar
    %Return the jj-th sample of parameters
    for dd=1:NStation+1
        par=cell2mat(params(dd,:));
        par=par(jj,:);
        PARR(dd,:)=par;%Concatenate the dd-th parameters 
    end
    %Calculation of the pseudovalues
    for dd=1:NStation
        J(dd,:)=((NStation).*PARR(NStation+1,:))-(NStation-1).*PARR(dd,:);%Concatenate the dd-th pseudovalue
    end
    %Calculation of the JackKnife estimator
    J_valor=mean(J);
    %Calculation of the standard error of the JackKnife estimator
    J_diff=zeros(NStation,4);%Initialize the difference matrix
    dd=0;
    for dd=1:NStation
        J_diff(dd,:)=(J(dd,:)-J_valor).^2;%Concatenate the dd-th difference
    end
    J_error(jj,:)=sqrt((1/(NStation*(NStation-1))).*sum(J_diff));%Concatenate the jj-th error values 
end
end