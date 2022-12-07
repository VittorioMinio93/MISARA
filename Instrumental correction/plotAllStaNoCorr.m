%%%%%%%%%%%%%Return all seismic traces%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pathname: the Input directory
% filez: the uncorrected file (.mat)
% file_name: the term of research file 
% text1: Window message
% comp_chan: the station system used (component/channel)
% comp: the name of the component/channel
% coord: the path of the coordinate file (.mat) 
% station, stations: station names
% nsta: number fo stations used
% sta_ref: station name of reference 
% ff: list of files with their path
% ff2: list of seismic traces directories
% mytrace: seismic trace structure
% dt: sampling interval
% fs: sampling rate
% nyq: nyquist number
% z: amplitude vector of the seismic trace
% nt: number of samples
% startTime, endTime: temporal limits of the traces 
% traceout, TraceOut: uncorrected seismic trace vector/matrix
% time, TimeOut: time vector/matrix

function [TraceOut,TimeOut,nsta,stations]=plotAllStaNoCorr(pathname,filez,coord,comp,comp_chan,text1)
%Create the term of research files
file_name=filez(end-18:end-4);
station=[];%Initialize the station names vector
%Load the station coordinates file
load(coord,'station')
%Error control
if isempty(coord) | isempty(station);set(text1, 'String','Error. No coordinate file directory or No correct coordinate file.');drawnow;return;end
stations=station;
nsta=length(stations);
TraceOut=cell(1,length(stations));%Initialize the seismic trace matrix
TimeOut=cell(1,length(stations));%Initialize the time matrix
%Search all seismic traces based on the station names
list=searchAllStaz(pathname,station,comp,file_name,comp_chan);
%Error control
if isempty(list);set(text1, 'String','Error. No files are found.');drawnow;return;end
nlist=length(list);%Number of traces selected
%Loop through the number of stations used
for kk=1:nlist
    %Load the kk-th seimic trace
    file1=strcat(pathname,list(kk).name);
    load(file1)
    z=mytrace.data;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate; nyq=0.5/dt; nt=mytrace.sampleCount;
    startTime=mytrace.startTime;endTime=mytrace.endTime;
    %Redefine the station name of reference
    switch comp_chan
        case 'Comp'
            sta_ref=mytrace.station;
        case 'Chan'
            sta_ref=[mytrace.station '.' mytrace.channel ];
    end
    idx=strcmp(stations,sta_ref);
    %Rescale and remove the best straight-fit line
    z=z-mean(z);
    z=detrend(z);
    %Create the time vector
    Time1=datetime(startTime,'ConvertFrom','datenum');
    Time2=datetime(endTime,'ConvertFrom','datenum');
    Time1.Format = 'dd-MM-yyyy HH:mm:ss.SSSSS';
    Time2.Format = 'dd-MM-yyyy HH:mm:ss.SSSSS';
    time=Time1+seconds(dt):seconds(dt):Time2;
    time=datenum(time');
    %Concatenate the results
    TraceOut{idx}=z;
    TimeOut{idx}=time;
end
end