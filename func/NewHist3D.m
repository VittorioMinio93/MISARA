%%%%%%%%%%%%%Calculation and plot the 3D Multivariate analysis%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% PAR, Dx: first parameter vector
% X: first parameter bins
% par, Dy: second parameter vector
% Y: second parameter bins
% PAR1, PAR2, PARstep: minimum, maximum and step first parameter value
% par1, par2, parstep: minimum, maximum and step second parameter value
% H: normalized bin counts matrix

function NewHist3D(PAR,par,PAR1,PAR2,PARstep,par1,par2,parstep)
%% Calculation of 3D the histograms
X=[PAR1:PARstep:PAR2]';
Y=[par1:parstep:par2]';
nx=length(X);ny=length(Y);%Get the length of parameters vectors
Dx=PAR;Dy=par;
idx=isnan(Dy);Dx(idx)=[];Dy(idx)=[];%Find and delete NaN values from the parameters vectors
if isempty(Dx); return;end%Stop the routine if the first parameter vector is empty
n = length(Dx) ;%Eventually, get the new length of the first parameter vector
H = zeros(ny,nx) ;%Initialize the normalized bin counts matrix
%loop through the parameters vectors 
for i = 1:n
    x = dsearchn(X,Dx(i)) ;%Search the index of the i-th value of the first parameter vector Dx in the first parameter bins X
    y = dsearchn(Y,Dy(i)) ;%Search the index of the i-th value of the second parameter vector Dy in the second parameter bins Y
    H(y,x) = H(y,x) + 1 ;%Count for the x, y indexes
end
H=H./sum(sum(H));%Normalize the bin counts so that the sum is 1
idx=isnan(H); H(idx)=0;%Find and replace NaN values from the normalized bin counts matrix 
%% Plotting of the 3D histograms
 surf(X,Y,H);%Plot the 3D histogram
view(2); colormap parula;
colorbar;
shading interp ;
end