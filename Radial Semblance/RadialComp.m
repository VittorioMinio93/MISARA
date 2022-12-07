%%%%%%%%%%%%%%%Calculation of the radial component%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Z,N,E:  the three components matrixs
% Xe,Ye,Ze: coordinates of the search grid point (UTM system)
% xs, ys, zs: station coordinates (UTM system)
% NStation: number of stations used
% M: number of samples
% deltax,deltay,deltaz: differences between stations and grid point coordinates
% inc: incidence (degrees)
% baz: backazimuth (degrees)
% L,Q,T: radial and transverse components
% sigmaL,sigmaQ,sigmaT: rms of the radial and transverse components 
% sigma: rms composition
% data: radial component matrix

function [data,sigma]=RadialComp(Xe,Ye,Ze,xs,ys,zs,Z,E,N,NStation,M)
%Loop through the number of stations used
for kk=1:NStation
    %Compute the distances between stations and grid point coordinates
    deltax=Xe-xs(kk);
    deltaz=zs(kk)-Ze;
    deltay=Ye-ys(kk);
    r=sqrt((Xe-xs(kk))^2+(Ye-ys(kk))^2);
    %Compute backazimuth and incidence
    inc=atand(r/deltaz);
    baz=90-atan2(deltay,deltax).*180./pi;
    k=find(baz<0); baz(k)=baz(k)+360;
    k2=find(inc<0); inc(k2)=inc(k2)+180;
    %Convert Z,N,E componets into L,Q,T components
    [L,Q,T]=signalrotation3D(Z(:,kk),E(:,kk),N(:,kk),baz,inc);
    %Calculation of the rms
    sigmaL=sqrt((1/M)*(sum(L.^2,1)));
    sigmaQ=sqrt((1/M)*(sum(Q.^2,1)));
    sigmaT=sqrt((1/M)*(sum(T.^2,1)));
    %Concatenate the results
    sigma(:,kk)=sqrt(sigmaL^2+sigmaQ^2+sigmaT^2);
    data(:,kk)=L;
end
end