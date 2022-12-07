%%%%%%%%%%%%% Beam Pattern algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xs, ys: coordinates in UTM system
% smin: the minimum slowness for the Beam pattern analysis (s/km)
% smax: the miximum slowness for the Beam pattern analysis (s/km)
% ns: the dimensions of the square grid slowness for the Beam pattern analysis
% nsta: number of stations used
% omega: frequency of the wavefield
% sx, sy: slowness grid vectors
% B: Beam Pattern matrix
function [B,sx,sy]=BeamPatternS(xs,ys,smin,smax,f,ns)
%Set the dimensions of the slowness grid
sx=linspace(smin,smax,ns);sy=linspace(smin,smax,ns);
nsta=length(xs);
ns=length(sx);
%%Calculate the station coordinates realtive to array center
xs=xs-mean(xs);ys=ys-mean(ys);
omega=2*pi*f;
%Loop through the slowness grid
for isx=1:ns
    for isy=1:ns
        %Calculation of the Beam Pattern function in every Grid point
        a=exp(-i*omega*([xs',ys']*[sx(isx);sy(isy)]));
        A(isx,isy,:)=a;
    end
end
%Normalize the results
B=(1/ns^2).*abs(sum(A,3)).^2;
end