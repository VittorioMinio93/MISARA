%%%%%%%%%%%%%Calculation and plot the 3D temporal histograms%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% time, Dx: time vector
% X: time bins
% par, Dy: parameter vector
% Y: parameter bins
% time1, time2, timestep: minimum, maximum and step time value
% par1, par2, parstep: minimum, maximum and step parameter value
% H: normalized bin counts matrix

function NDistrHist3D(time,par,time1,time2,timestep,par1,par2,parstep)
%% Calculation of 3D the temporal histograms
X=[time1:datenum(timestep):time2]';
Y=[par1:parstep:par2]';
nx=length(X);ny=length(Y);%Get the length of the time and parameter vectors
Dx=time;Dy=par;
idx=isnan(Dy);Dx(idx)=[];Dy(idx)=[];%Find and delete NaN values from the time and the parameter vectors
if isempty(Dx); return;end%Stop the routine if the time vector is empty
n = length(Dx) ;%Eventually, get the new length of the time vector
H = zeros(ny,nx) ;%Initialize the normalized bin counts matrix
%loop through the time and parameter vectors 
for i = 1:n
    x = dsearchn(X,Dx(i)) ;%Search the index of the i-th value of the time vector Dx in the time bins X
    y = dsearchn(Y,Dy(i)) ;%Search the index of the i-th value of the parameter vector Dy in the parameter bins Y
    H(y,x) = H(y,x) + 1 ;%Count for the x, y indexes
end
X=datetime(X,'ConvertFrom','datenum');%Convert the time bins in datetime format
H=H./sum(H);%Normalize the bin counts so that the sum is 1
idx=isnan(H); H(idx)=0;%Find and replace NaN values from the normalized bin counts matrix 
%% Plotting of the 3D temporal histograms
 surf(X,Y,H);%Plot the 3D histogram
view(2); colormap('parula'); 
colorbar
shading ('interp') 
end