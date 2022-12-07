%%%%%%%%%%%%%Save the corrected seismic traces%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pathname: the Input directory
% save2, pathout: the Output directories
% filez: the uncorrected file (.mat)
% file_name: the term of research file 
% param: the instrument response correction parameters file (.mat)
% text1: Window message
% comp_chan: the station system used (component/channel)
% component: the name of the component/channel
% coord: the path of the coordinate file (.mat) 
% station, stations: station names
% Data: the table of the instrument response correction parameters 
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
% pols: instrument Poles from factory calibration worksheet (rad/s) 
% zers: instrument Zeroes from factory calibration worksheet (rad/s) 
% bmv : coefficient for conversion from bits to Volts (V/counts) 
% vvel: coefficient for conversion from Volts to m/s (V s/m)
% k = calibration coefficient (rad/s)
% traceout: corrected seismic velocity trace

function savenewfile(pathname,filez,param,save2,text1,comp_chan,component,coord)
station=[];%Initialize the station names vector
try
    %Create the term of research files
    file_name=filez(end-18:end-4);
    pathout=[];%Initialize the Output vector
    %Load the station coordinates file
    load(coord)
    %Error control
    if isempty(coord) | isempty(station);set(text1, 'String','Error. No coordinate file directory or No correct coordinate file.');drawnow;return;end
    %Load the table of the instrument response correction parameters
    load(param)
    stations=Data.stations;
    %Search all seismic traces based on the station names
    list=searchAllStaz(pathname,station,component,file_name,comp_chan);
    %Error control
    if isempty(list);set(text1, 'String','Error. No files are found.');drawnow;return;end
    nlist=length(list);%Number of traces selected
    %Loop through the number of stations used
    for kk=1:nlist
        %Load the kk-th seimic trace
        file1=strcat(pathname,list(kk).name);
        load(file1)
        z=mytrace.data;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate; nyq=0.5/dt; nt=mytrace.sampleCount;
        startTime=mytrace.startTime;endTime=mytrace.endTime;Time1=datetime(startTime,'ConvertFrom','datenum');
        %Redefine the station name of reference
        switch comp_chan
            case 'Comp'
                sta_ref=mytrace.station;
            case 'Chan'
                sta_ref=[mytrace.station '.' mytrace.channel ];
        end
        %Select the analysis parameters
        idx=strcmp(stations,sta_ref);
        pols=Data.pols{idx};
        zers=Data.zers{idx};
        bmw=Data.bmw(idx);
        k=Data.k(idx);
        vvel=Data.vvel(idx);
        jday=num2str(day(Time1,'dayofyear'));
        %Rescale and remove the best straight-fit line
        z=z-mean(z);
        z=detrend(z);
        %Create the Output directory
        pathout=strcat(save2,'/COR_SIGN/',num2str(jday),'/');if isempty(pathout);mkdir(pathout);end
        %Apply the instrument response correction
        traceout = instrument_resp_remove(nyq,pols,zers,k,bmw,vvel,z);
        mytrace.data=traceout;
        %Save the corrected seismic trace
        save([pathout list(kk).name],'mytrace')
    end
    %Command message
    set(text1, 'String', sprintf('Inst_Corr.mat-Observed data are saved.'));
    drawnow
catch
    %Error assessment
    if isempty(z);set(text1, 'String','Error. No correct signals data.');drawnow;return;end
end
end