%%%%%%%%%%%%%%%%%%%%%%Radial Semblance algorithm%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Z,N,E: the three components matrixs
% ZZ,NN,EE: the cutted three components matrixs
% M: number of samples
% data: the radial component maxtrix
% sigma: the rms of the three components data 
% Xq,Yq,Zq: 3D search grid Nx*Ny*Nz (UTM system)
% xq,yq,zq: coordinate vectors of the search grid (UTM system)
% loc: 3D Radial Semblance grid Nx*Ny*Nz
% xs, ys, zs: station coordinates (UTM system)
% NStation: number of stations used
% out: distances from stations positions (km)
% v: the velocity value (km/s)
% Nx,Ny,Nz: dimensions of the search Grid
% fs: sampling rate
% trig: picking time reference (samples)
% staidx: station index of reference
% wl: analysis window (samples)
% valor: maximum Radial Semblance value
% ix1,ix2,ix3: indexes of the maximum Radial Semblance value

function [valor,ix1,ix2,ix3,loc]=RadialS(Z,E,N,Xq,Yq,Zq,xq,yq,zq,NStation,v,Nx,Ny,Nz,xs,ys,zs,fs,trig,staidx,w1)
loc = zeros(Ny,Nx,Nz);%Initialize the Radial Semblance matrix
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
                %Initialize the cutted three components matrixs
                ZZ=zeros(w1,NStation);
                EE=zeros(w1,NStation);
                NN=zeros(w1,NStation);
                %Cut the three components of the signal of reference
                ZZ(:,staidx)=Z(t1:t2,staidx);
                EE(:,staidx)=E(t1:t2,staidx);
                NN(:,staidx)=N(t1:t2,staidx);
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
                        %Time shifting of the three components
                        ZZ(:,ii)=Z(t1+DDelay:t2+DDelay,ii);
                        EE(:,ii)=E(t1+DDelay:t2+DDelay,ii);
                        NN(:,ii)=N(t1+DDelay:t2+DDelay,ii);
                    end
                end
                M=size(ZZ,1);
                % Calculation of the radial component of the signals
                [data,sigma]=RadialComp(xq(nnx),yq(nny),zq(nnz),xs,ys,zs,ZZ,EE,NN,NStation,M);
                %Normalize the radial component of the signals
                data=data./sigma;
                %Application of the Radial Semblance algorithm
                loc(nny,nnx,nnz)=sum(sum(data,2).^2+NStation*sum(data.^2,2))/(2*M*NStation^2);
            end
        end
    end
end

%Points of maximum probability
[valor,ix]=max(loc(:),[],'omitnan');
[ix1, ix2,ix3] = ind2sub(size(loc),ix);
end
