%%%%%%%%%%%%%Jurkevics algorithm%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% az: azimuth vectors
% incid: incidence vectors
% rect: rectilinearity vectors
% sis: three components matrix 
% wl, ws: windows analysis (s)
% wls, wss: windows analysis (samples) 
% f1, f2: frequency analysis (Hz)
% fs: sampling rate
% dt: sampling interval
% nyq: nyquist number
% temp: covariance matrix
% U: eigenvectors
% Dt: eigenvalues
% z,n,e: the three direction cosines of the eigenvector corrisponding with the largest eigenvalue
% LL: sum of the auto-variances  
% TL: time samples

function  [az,incid,rect,LL,TL] =  polarization(sis,fs,wl,ws)
%Set the analysis parameters
dt=1/fs;    
nyq=fs/2;		
ns=length(sis);%The row dimension of three components matrix 
t=0:dt:(ns-1)*dt;% time samples
wls=round(wl/dt);
wss=round(ws/dt);
nw=1+fix((ns-wls)/wss);%number of analysis window
%Move the analysis window through the seismic traces
for k=1:nw
    i1=1+(k-1)*wss;i2=i1+wls-1;%Limits of the analysis window
    %Calculation of covariance matrix 
    temp=cov(sis(i1:i2,:));
    %Control the calculation through the covariance values
    if (length(find(temp==0)) == 0 )
        %Eigenvalues and eigenvectors
        [U Dt]=eig(temp);
        %Creation of the diagonal matrix
        D=diag(Dt);
        %Estimation of the main axes of the polarization ellipsoid
        z=U(1,3);
        n = U(2,3)*sign(z);
        e = U(3,3)*sign(z);
        r = n/e;
        %Calcultion of the azimuth from the eigenvectors corrisponding with the largest eigenvalue
        %Set the excepts
        if ( r == 0 ) 
            if z > 0 
                az(k) = 0;
            else
                az(k) = 180;
            end
        else
            if ( r == inf )
                if z > 0 
                    az(k) = 90;
                else
                    az(k) = 270;
                end
            else
                az(k)= 90-(180/pi)*atan(n./e);
                %Correction based on the sign
                if ( (n < 0 ) && (e < 0) ) || ( ( n > 0 ) && ( e < 0 ) )
                    az(k) = az(k) + 180;
                end
            end
        end
        %Calcultion of the incidence from the eigenvectors corrisponding with the largest eigenvalue
        chor=sqrt(U(2,3).^2+U(3,3).^2);
        inrad=atan(chor./U(1,3));
        incid(k)=(180/pi).*(inrad*sign(inrad));
        %Calcultion of the rectilinearity from eigenvalues
        L1 = D(1);
        L2 = D(2);
        L3 = D(3);
        NL =  (1-L2/L3)^2 + (1-L1/L3)^2 + ( (L2/L3) - (L1/L3) )^2;
        DL = 2*((1 + L2/L3 - L1/L3))^2;
        rect(k) = sqrt(NL/DL);
    else
        az(k) = NaN;
        rect(k) = NaN;
        incid(k) = NaN;
    end
    LL(k)=sum(diag(temp));
    TL(k)=((i1+i2)/2).*dt;
end
end