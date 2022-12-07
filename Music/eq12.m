%%%%%%%%%%%%%Evalutate power spectrum over a grid of slownesses%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N: number of stations
% freq: the frequency analysis
% E: eigenvectors selected matrix
% L: projection length
% x: station coordinates matrix
% slowx0, slowy0: the limits of the grid slowness (s/km)
% ds: the size step of the grid slowness (s/km) 
% ns: dimensions of slowness grid
% flagsig: method to determine ranges of eigenvectors to use
% sx_shift, sy_shift: the initial slowness
% Nsig: the number of sources computed 
% Nsig2,Nsig2: number of signal sources/ signal noise
% UE: the directional function
% power: the power spectrum matrix

function power=eq12(N,freq,E,L,x,row1,slowx0,slowy0,ds,ns,flagsig,sx_shift,sy_shift,Nsig,Nsig1,Nsig2)
pi=4*atan(1);
power=zeros(ns,ns);%Initialize the power spectrum matrix
%Loop over slowness grid
slowy = slowy0-sy_shift;
for i=1:ns
    ky = slowy*freq*2.*pi;
    slowx = slowx0-sx_shift;
    for ii=1:ns
        kx = slowx*freq*2.*pi;
        %Compute the directional function and the power spectrum
        UE=UdotE(N,kx,ky,E,L,x,row1,flagsig,Nsig,Nsig1,Nsig2);
        power(ii,i) = real(UE);
        slowx = slowx + ds;
    end
    slowy = slowy + ds; 
end
%Rescale power spectrum
if (flagsig=='H') | (flagsig=='N')
    for i=1:ns
        for j=1:ns
            power(i,j) = 1./power(i,j);
        end
    end   
end
end