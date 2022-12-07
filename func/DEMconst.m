%%%%%%%%%%%%% Construct the Search grid and the DEM grid%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dem_file: path of the DEM file variables
% correction_topo: the topographic correction (yes o no)
% xmin: the minimum X-limit of the search grid  (km)
% xmax: the maximum X-limit of the search grid  (km)
% ymin: the minimum Y-limit of the search grid  (km)
% ymax: the maximum Y-limit of the search grid  (km)
% zmin: the minimum Z-limit of the search grid  (km)
% zmax: the maximum Z-limit of the search grid  (km)
% xq, yq, zq: coordinate vectors of the search grid (UTM system)
% step_sem: the search grid step size (km)
% X,Y,A: georeferenced DEM grid NxM (UTM system)
% xLimits, yLimits: georeferenced limits of the DEM grid in m
% xdem_UTM, ydem_UTM: coordinate vectors of the DEM grid (km) 
% Xq, Yq, Zq: 3D search grid Nx*Ny*Nz (UTM system)
%Perfil_A: Interpolated elevation grid (km)

function [Xq,Yq,Zq,Nx,Ny,Nz,xq,yq,zq,A,X,Y]=DEMconst(dem_file,correction_topo,xmin,xmax,ymin,ymax,zmin,zmax,step_sem)
%% Grid Construction
load(dem_file)%load A, xLimits, yLimits variables
[N,M]= size(A);%Get the dimensions of the elevetion grid NxM
ydem_UTM= linspace(yLimits(2), yLimits(1), N)/1000;
xdem_UTM= linspace(xLimits(1), xLimits(2), M)/1000;
xq=xmin:step_sem:xmax;
yq=ymin:step_sem:ymax;   
zq=zmin:step_sem:zmax;
Nx=length(xq); Ny=length(yq); Nz=length(zq);%Get the dimensions of the search grid
[Xq,Yq]=meshgrid(xq,yq);%Get the 2D search grid Nx*Ny
[X,Y] = meshgrid(xdem_UTM, ydem_UTM);%Get georeferenced DEM grid NxM
Perfil_A = interp2(X,Y,A,Xq,Yq);%Interp the DEM grid to the search grid
[Xq,Yq,Zq]=meshgrid(xq,yq,zq);%Get the 3D search grid Nx*Ny*Nz
%% Apply the topographic correction to the search grid 
if strcmp(correction_topo,'yes')
    %loop  through every surface of the Zq 3D search grid (TT)
        for vv=1:Nz
             TT=Zq(:,:,vv);
             cond=TT>Perfil_A;%Find the grid position located over the topographic surface
             TT(cond)=NaN;
             Zq(:,:,vv)=TT;%Assign a NaN value to all grid position over the topographic surface
        end
end
end