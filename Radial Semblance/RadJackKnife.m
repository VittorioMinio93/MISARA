%%%%%%%%%%%JackKnife method for the Radial Semblance algorithm%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Z,N,E: the three components matrixs
% Z_1,N_1,E_1: the "leave one out" three components matrixs
% ZZ,NN,EE: the cutted three components matrixs
% M: number of samples
% data: the radial component maxtrix
% sigma: the rms of the three components data 
% Xq,Yq,Zq: 3D search grid Nx*Ny*Nz (UTM system)
% xq,yq,zq: coordinate vectors of the search grid (UTM system)
% loc: 3D Radial Semblance grid Nx*Ny*Nz
% xs, ys, zs: station coordinates (UTM system)
% xs_2, ys_2, zs_2: "leave one out" station coordinates (UTM system)
% NStation: number of stations used
% out,out_2: distances from stations positions (km)
% v: the velocity value (km/s)
% Nx,Ny,Nz: dimensions of the search Grid
% fs: sampling rate
% trig: picking time reference (samples)
% staidx: station index of reference
% wl: analysis window (samples)
% valor: maximum Radial Semblance value
% ix1,ix2,ix3: indexes of the maximum Radial Semblance value
% ixx1,ixx2,ixx3: indexes of the maximum Radial Semblance value
% Coord: source locations  computed by omitting one station  at a time (UTM system)
% J: pseudovalues
% J_valor: JackKnife estimator of the parameters computed by omitting one station 
% J_diff: squared difference between the pseudovalues and the JackKnife estimator 
% J_error: standard error of the JackKnife estimator (km)

function [J_error]=RadJackKnife(Z,E,N,Xq,Yq,Zq,xq,yq,zq,NStation,v,Nx,Ny,Nz,xs,ys,zs,fs,trig,staidx,w1,ixx1,ixx2,ixx3)
out(NStation) = struct('distance',zeros(Ny,Nx,Nz));%Initialize the distance matrix
%Set the limits of the analysis window
t1=trig-w1+1;
t2=trig;
%Compute the distance between the search grid points and the stations positions 
for pp=1:NStation
    out(pp).distance=(sqrt( (xs(pp)-Xq).^2 +  (ys(pp)-Yq).^2 +  (zs(pp)-Zq).^2 ));
end
%Loop through the number of station used
for gg=NStation:-1:1
    loc = zeros(Ny,Nx,Nz);%Initialize the Radial Semblance matrix
    out_2(NStation-1) = struct('distance',zeros(Ny,Nx,Nz));%Initialize the distance matrix omitting the gg-th station
    nny=0;
    nnx=0;
    nnz=0;
    %Initialize the distance matrix omitting the gg-th station
    Z_1=Z;
    Z_1(:,gg)=[];
    E_1=E;
    E_1(:,gg)=[];
    N_1=N;
    N_1(:,gg)=[];
    xs_2=xs;
    xs_2(gg)=[];
    ys_2=ys;
    ys_2(gg)=[];
    zs_2=zs;
    zs_2(gg)=[];
    pp=0; 
    %Compute the distance between the search grid points and the stations positions 
    for pp=1:NStation-1
        out_2(pp).distance=(sqrt( (xs_2(pp)-Xq).^2 +  (ys_2(pp)-Yq).^2 +  (zs_2(pp)-Zq).^2 ));
    end
    %Loop through the 3D search grid
    for nny = 1:Ny
        for nnx = 1:Nx
            for nnz = 1:Nz
                %Return NaN values for  NaN distances
                if isnan(out(1).distance(nny,nnx,nnz))
                    loc(nny,nnx,nnz)=NaN;
                else
                    %Control the routines on the basis of the station index of reference
                    if gg>staidx
                        %Initialize the cutted three components matrixs
                        ZZ=zeros(w1,NStation-1);
                        EE=zeros(w1,NStation-1);
                        NN=zeros(w1,NStation-1);
                        %Cut the three components of the signal of reference
                        ZZ(:,staidx)=Z_1(t1:t2,staidx);
                        EE(:,staidx)=E_1(t1:t2,staidx);
                        NN(:,staidx)=N_1(t1:t2,staidx);
                        %Loop through the new number of stations used
                        for ii = 1:NStation-1
                            if ii~=staidx
                                %Get the distance of the ii-th station from the grid point
                                r=out_2(ii).distance(nny,nnx,nnz);
                                %Get the distance of the station of reference from the grid point
                                rmain=out(staidx).distance(nny,nnx,nnz);
                                %Calculation of delay time
                                NDelay=round((r/v)*fs);
                                MDelay=round((rmain/v)*fs);
                                DDelay=NDelay-MDelay;
                                %Time shifting of the three components
                                ZZ(:,ii)=Z_1(t1+DDelay:t2+DDelay,ii);
                                EE(:,ii)=E_1(t1+DDelay:t2+DDelay,ii);
                                NN(:,ii)=N_1(t1+DDelay:t2+DDelay,ii);
                            end
                        end
                        M=size(ZZ,1);
                        %Calculation of the radial component of the signals
                        [data,sigma]=RadialComp(xq(nnx),yq(nny),zq(nnz),xs_2,ys_2,zs_2,ZZ,EE,NN,NStation-1,M);
                        %Normalize the radial component of the signals
                        data=data./sigma;
                        %Application of the Radial Semblance algorithm
                        loc(nny,nnx,nnz)=sum(sum(data,2).^2+(NStation-1)*sum(data.^2,2))/(2*M*(NStation-1)^2);
                    elseif gg==staidx
                        %Initialize the cutted three components matrixs
                        ZZ=zeros(w1,NStation-1);
                        EE=zeros(w1,NStation-1);
                        NN=zeros(w1,NStation-1);
                        for ii = 1:NStation-1
                            %Get the distance of the ii-th station from the grid point
                            r=out_2(ii).distance(nny,nnx,nnz);
                            %Get the distance of the station of reference from the grid point
                            rmain=out(staidx).distance(nny,nnx,nnz);
                            %Calculation of delay time
                            NDelay=round((r/v)*fs);
                            MDelay=round((rmain/v)*fs);
                            DDelay=NDelay-MDelay;
                            %Time shifting of the three components
                            ZZ(:,ii)=Z_1(t1+DDelay:t2+DDelay,ii);
                            EE(:,ii)=E_1(t1+DDelay:t2+DDelay,ii);
                            NN(:,ii)=N_1(t1+DDelay:t2+DDelay,ii);   
                        end
                        M=size(ZZ,1);
                        %Calculation of the radial component of the signals
                        [data,sigma]=RadialComp(xq(nnx),yq(nny),zq(nnz),xs_2,ys_2,zs_2,ZZ,EE,NN,NStation-1,M);
                        %Normalize the radial component of the signals
                        data=data./sigma;
                        %Application of the Radial Semblance algorithm
                        loc(nny,nnx,nnz)=sum(sum(data,2).^2+(NStation-1)*sum(data.^2,2))/(2*M*(NStation-1)^2);
                    elseif gg<staidx  
                        %Initialize the cutted three components matrixs
                        ZZ=zeros(w1,NStation-1);
                        EE=zeros(w1,NStation-1);
                        NN=zeros(w1,NStation-1);
                        %Cut the three components of the signal of reference
                        ZZ(:,staidx-1)=Z_1(t1:t2,staidx-1);
                        EE(:,staidx-1)=E_1(t1:t2,staidx-1);
                        NN(:,staidx-1)=N_1(t1:t2,staidx-1);
                        for ii = 1:NStation-1
                            if ii~=staidx-1
                                %Get the distance of the ii-th station from the grid point
                                r=out_2(ii).distance(nny,nnx,nnz);
                                %Get the distance of the station of reference from the grid point
                                rmain=out(staidx).distance(nny,nnx,nnz);
                                %Calculation of delay time
                                NDelay=round((r/v)*fs);
                                MDelay=round((rmain/v)*fs);
                                DDelay=NDelay-MDelay;
                                %Time shifting of the three components
                                ZZ(:,ii)=Z_1(t1+DDelay:t2+DDelay,ii);
                                EE(:,ii)=E_1(t1+DDelay:t2+DDelay,ii);
                                NN(:,ii)=N_1(t1+DDelay:t2+DDelay,ii);
                            end
                        end
                        M=size(ZZ,1);
                        %Calculation of the radial component of the signals
                        [data,sigma]=RadialComp(xq(nnx),yq(nny),zq(nnz),xs_2,ys_2,zs_2,ZZ,EE,NN,NStation-1,M);
                        %Normalize the radial component of the signals
                        data=data./sigma;
                        %Application of the Radial Semblance algorithm
                        loc(nny,nnx,nnz)=sum(sum(data,2).^2+(NStation-1)*sum(data.^2,2))/(2*M*(NStation-1)^2);
                    end
                end
            end
        end
    end
    
    %Points of max probability
    [valor,ix]=max(loc(:),[],'omitnan');
    [ix1, ix2,ix3] = ind2sub(size(loc),ix);
    %Calculating of the source coordinates omitting the gg-th station
    Coord_x=xq(ix2);
    Coord_y=yq(ix1);
    Coord_z=zq(ix3);
    Coord_xyz=[Coord_x,Coord_y,Coord_z];
    %Concatenate the JackKnife results
    Coord(gg,:)=Coord_xyz;
    %Clear the localization variables 
    clear loc ys_2 xs_2 zs_2 Z_1 E_1 N_1
end
%Calculating of the source coordinates for all stations
Coord_x=xq(ixx2);
Coord_y=yq(ixx1);
Coord_z=zq(ixx3);
Coord_xyz=[Coord_x,Coord_y,Coord_z];
Coord(NStation+1,:)=Coord_xyz;
%Application of the  Jacknife method on the source locations
J=zeros(NStation,3);%Initialize the pseudovalue vector
%Calculation of the pseudovalues
for dd=1:NStation
    J(dd,:)=((NStation).*Coord(NStation+1,:))-(NStation-1).*Coord(dd,:);
end
%Calculation of the JackKnife estimator
J_valor=mean(J);
%Calculation of the standard error of the JackKnife estimator
J_diff=zeros(NStation,3);%Initialize the difference matrix
dd=0;
for dd=1:NStation
    J_diff(dd,:)=(J(dd,:)-J_valor).^2;%Concatenate the dd-th difference
end
J_error=sqrt((1/(NStation*(NStation-1))).*sum(J_diff));
end