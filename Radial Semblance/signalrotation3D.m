%%%%%%%%%%%%%%%Three dimensional rotation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Z,N,E: the three components matrixs
% inc: incidence (degrees)
% baz: backazimuth (degrees)
% L,Q,T: radial and transverse components 

function [L,Q,T]=signalrotation3D(Z,E,N,ba,inc)
%Convert Z,N,E componets into L,Q,T components
L=Z*cosd(inc)-N*sind(inc)*cosd(ba)-E*sind(inc)*sind(ba);
Q=Z*sind(inc)+N*cosd(inc)*cosd(ba)+E*cosd(inc)*sind(ba);
T = N*sind(ba)-E*cosd(ba); 
end