%%%%%%%%%%%%% Calculation of the kinematic attributes%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X: station coordinates matrix
% az: azimuth vector
% rayp: ray paramater vector
% m: slowness matrix
% sx: horizontal slowness vector
% sy: horizontal slowness vector
% nsta: number of station used
% delay: delay time 

function [sx sy az rayp]=FitPlane(X,delay)
nsta=length(X);
ip=0;
%Loop through all independent pairs of stations
for is1=1:nsta-1
    for is2=is1+1:nsta
        %Calculation of the differences between the station coordinates
        ip=ip+1;
        G(ip,1)=X(is2,1)-X(is1,1);%Concatenate the ip-th x-coordinate
        G(ip,2)=X(is2,2)-X(is1,2);%Concatenate the ip-th y-coordinate
    end
end
[np nt]=size(delay);%dimensions of delay time matrix
sx=zeros(nt,1);%Initialize the x-component vector of the slowness
sy=zeros(nt,1);%Initialize the y-component vector of the slowness
%Loop through the number of samples
for k=1:nt
    %Calculation of the slowness matrix through the generalized inverse
    m=(inv(G'*G))*G'*(delay(:,k));
    %Concatenate the horizontal components of the slowness
    sx(k)=m(1);
    sy(k)=m(2);
end
%Calculation of the azimuth from the horizontal slowness
az=90-atan2(sy,sx).*180./pi;
k=find(az<0); az(k)=az(k)+360;
%Calculation of the ray parameter from the horizontal slowness
rayp=sqrt(sx.^2+sy.^2);
end