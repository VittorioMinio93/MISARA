%%%%%%%%%%%%%%%%%%%%%%Semblance algorithm%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Z: the signals matrix
% ZZ,data: the cutted signals matrix
% Xq,Yq,Zq: 3D search grid Nx*Ny*Nz (UTM system)
% loc: 3D Semblance grid Nx*Ny*Nz
% xs, ys, zs: station coordinates (UTM system)
% NStation: number of stations used
% out: distances from stations positions (km)
% v: the velocity value (km/s)
% Nx,Ny,Nz: dimensions of the search Grid
% fs:sampling rate
% trig:picking time reference (samples)
% staidx: station index of reference
% wl: analysis window (samples)
% K: the frequency-dependent absorption coefficient for the seismic amplitude decays with distance
% n: the exponent value 
% valor: maximum Semblance value
% ix1,ix2,ix3: indexes of the maximum Semblance value

function [valor,ix1,ix2,ix3,loc]=Semblance(Z,Xq,Yq,Zq,NStation,v,Nx,Ny,Nz,xs,ys,zs,fs,trig,staidx,w1,K,n)
loc = zeros(Ny,Nx,Nz); %Initialize the Semblance matrix
out(NStation) = struct('distance',zeros(Ny,Nx,Nz));%Initialize the distance matrix
%Set the limits of the analysis window
t1=trig-w1+1;
t2=trig;
%Compute the distance between the search grid points and the stations positions 
for pp=1:NStation
    out(pp).distance=(sqrt( (xs(pp)-Xq).^2 +  (ys(pp)-Yq).^2 +  (zs(pp)-Zq).^2 ));
end
%Loop through the 3D search grid
for nny = 1:Ny
    for nnx = 1:Nx
        for nnz = 1:Nz
            %Return NaN values for  NaN distances
            if isnan(out(1).distance(nny,nnx,nnz))
                loc(nny,nnx,nnz)=NaN;
            else
                ZZ=zeros(w1,NStation);%Initialize the cutted signals matrix
                %Cut the signal of reference
                ZZ(:,staidx)=Z(t1:t2,staidx);
                %Loop through the number of stations used
                for ii = 1:NStation
                    if ii~=staidx
                        %Get the distance of the ii-th station from the grid point
                        r=out(ii).distance(nny,nnx,nnz);
                        %Get the distance of the station of reference from the grid point
                        rmain=out(staidx).distance(nny,nnx,nnz);
                        %Calculation of delay time
                        NDelay=round((r/v)*fs);
                        MDelay=round((rmain/v)*fs);
                        DDelay=NDelay-MDelay;
                        %Time shifting of the cutted signals
                        ZZ(:,ii)=Z(t1+DDelay:t2+DDelay,ii);
                        %Apply the amplitude decays with the distance r
                        ZZ(:,ii)=ZZ(:,ii)*(abs(r)^(n))*exp(K*abs(r));
                    end
                end
                %Apply the amplitude decays with the distance rmain
                ZZ(:,staidx)=ZZ(:,staidx)*(abs(rmain)^(n))*exp(K*abs(rmain));
                data=ZZ;
                %Application of the Semblance algorithm
                loc(nny,nnx,nnz)=(sum(sum(data,2).^2))/(NStation*sum(sum(data.^2,2)));
            end
        end
    end
end

%Points of maximum probability
[valor,ix]=max(loc(:),[],'omitnan');
[ix1, ix2,ix3] = ind2sub(size(loc),ix);
end
