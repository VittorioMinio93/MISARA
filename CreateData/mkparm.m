%%%%%%%%%%%%% Save instrumental response correction structures%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pathout: the Output folder 
% pathpar: the Output Instrument response directory
% net: network of reference
% K, k: calibration coefficient (rad/s)
% poles, pols: poles of the transfer function
% zeros, xers: zeros of the transfer function
% Vvel, vvel: coefficient for conversion from Volts to m/s (V s/m)
% Bmw, bmw: coefficient for conversion from counts to Volts (V/counts)
% sta: the name of the reference station 
% loc: localization of reference
% Data: table of the instrumental response correction parameters

function mkparm(pathout,net,K,poles,zeros,Vvel,Bmw,sta,loc) 
Data=[];%Initialize the table of the instrumental response correction parameters
%Create the Output directory
pathpar=strcat(pathout,'/',net,'/Instrument Response/');if ~exist(pathpar, 'dir');mkdir(pathpar);end
%Search the instrumental response file
p=dir(strcat(pathpar,'correction_parameters.mat'));
%Control and eventually set the instrumental response correction parameters
if isempty(K);K=NaN;end
if isempty(poles);poles=NaN;end
if isempty(zeros);zeros=NaN;end
if isempty(Vvel);Vvel=NaN;end
if isempty(Bmw);Bmw=NaN;end
staz=[sta loc];
if length(p)==0
    %Create the new instrumental response vectors 
   k=K;pols={poles};zers={zeros};vvel=Vvel;bmw=Bmw;stations={staz};
else
    %Load and Upload the new instrumental response values 
    load(strcat(pathpar,'correction_parameters.mat'))
    k=Data.k;bmw=Data.bmw;vvel=Data.vvel;
    pols=Data.pols;zers=Data.zers;stations=Data.stations;
    staidx=strcmp(staz,stations);
    k(staidx)=[];bmw(staidx)=[];vvel(staidx)=[];pols(staidx)=[];zers(staidx)=[];stations(staidx)=[];
     k=[k;K];pols=[pols;{poles}];zers=[zers;{zeros}];vvel=[vvel;Vvel];bmw=[bmw;Bmw];stations=[stations;{staz}];
end
%Save the table of the instrumental response correction parameters
Data=table(stations,pols,zers,bmw,k,vvel);Data.Properties.VariableNames={'stations','pols','zers','bmw','k','vvel'};
save([pathpar 'correction_parameters.mat'],'Data')
end