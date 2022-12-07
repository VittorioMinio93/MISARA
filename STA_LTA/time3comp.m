%%%%%%%%%%%%% Save the waveforms detected (component system) %%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sta_ref: the name of the reference station
% comp, component: the name of the components
% Timesel: trigger time (datenum format)
% Name: the path of the file of reference 
% path:the path of the stations coordinate file (.mat) 
% mainpath: the Input folder
% fs: sampling rate
% dt: sampling interval
% nt: number of samples
% wlong: the pre-trigger and post-trigger window (s)
% pathout: the Output folder
% outputformat:  the format (.mat/.txt) of the output files 
% text1:  Window message 
% ii: ii-th waveform
% bb: Total number of waveforms
% station: station names
% z,z3: amplitude vector of the seismic trace
% Z: one station waveform extracted
% mytrace: seismic trace structure
% TP_sampl: trigger time (samples)
% juliand, juliandnew: julian day

function time3comp(Timesel,Name,path,mainpath,fs,wlong,pathout,outputformat,text1,ii,bb)
try
    Z=[];%Initialize the waveform extracted vector 
    station=[];%Initialize the station names vector
    Timesel=datetime(Timesel,'ConvertFrom','datenum');%Convert in datetime format
    %Set the component names
    components={'Z','N','E'};
    %Load the station coordinates file 
    load(path,'station');
    %Error control
    if isempty(station);set(text1, 'String','Error. No coordinate file directory or No correct coordinate file.');drawnow;return;end
    %Loop through the three components
    for jj=1:length(components)
        comp=components{jj};%Select the jj-th component
        %Calculation of the trigger time in samples
        TP_sampl=(minute(Timesel)*60+round(second(Timesel)))*round(fs);
        %Compute the juliand day
        if day(Timesel,'dayofyear')<10
            juliand=['00' num2str(day(Timesel,'dayofyear'))];
        elseif day(Timesel,'dayofyear')>=10 & day(Timesel,'dayofyear')<100
            juliand=['0' num2str(day(Timesel,'dayofyear'))]; 
        elseif day(Timesel,'dayofyear')>=100 
            juliand=num2str(day(Timesel,'dayofyear')); 
        end
        %Loop through the number of stations used
        for kk=1:length(station)
            %Set the station of reference
            sta_ref=station{kk};
            %Create the name of the signal file 
            filesac=strcat(mainpath,'/',juliand,'/',sta_ref,'.BH',comp,Name(end-19:end));
            %Control the directory of the signal file
            p=dir(filesac);
            if length(p)~=0
                %Load the trace
                load(filesac); 
                z=mytrace.data;nt=mytrace.sampleCount;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate;
                %Rescale and remove the best straight-fit line
                z=z-mean(z);
                z=detrend(z);
                w1=round(wlong*fs);
                Z=[];%Initialize the waveform extracted vector 
                %Extracting the waveforms on the basis of the trigger time  and signal files available
                if TP_sampl-w1>0 & TP_sampl+w1<=length(z)
                    %Extraction of the waveform
                    Z=z(TP_sampl-w1:TP_sampl+w1);
                elseif TP_sampl-w1<0 & TP_sampl+w1<=length(z)
                    %Search the new time of reference
                    newtime=Timesel-1/24;
                    newtime2=datestr(newtime,'yyyymmdd-HHMMSS');
                    newtime2=newtime2(end-5:end-4);
                    %Compute the new juliand day
                    if day(newtime,'dayofyear')<10
                        juliandnew=['00' num2str(day(newtime,'dayofyear'))];
                    elseif day(newtime,'dayofyear')>=10 & day(newtime,'dayofyear')<100
                        juliandnew=['0' num2str(day(newtime,'dayofyear'))]; 
                    elseif day(newtime,'dayofyear')>=100 
                        juliandnew=num2str(day(newtime,'dayofyear')); 
                    end
                    %Create the new name of the signal file
                    newfile=strcat(mainpath,'/',juliandnew,'/',sta_ref,'.BH',comp,'.',num2str(year(newtime)),'.',juliandnew,'.',newtime2,'0000.sac');
                    %Control the directory of the new signal file 
                    p3=dir(newfile);
                    if length(p3)~=0
                        %Load the trace
                        load(newfile); 
                        z3=mytrace.data;nt=mytrace.sampleCount;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate;
                        %Rescale and remove the best straight-fit line
                        z3=z3-mean(z3);
                        z3=detrend(z3);
                        %Extraction of the waveform
                        Z1=z(1:TP_sampl+w1);
                        Z2=z3(end-(w1-TP_sampl):end);
                        Z=[Z2;Z1];
                    else
                        %Warning message
                        set(text1, 'String', 'Warning. One or more signals are no available.');
                        drawnow
                    end
                elseif TP_sampl-w1>0 & TP_sampl+w1>length(z)
                    %Search the new time of reference
                    newtime=Timesel+1/24;
                    newtime2=datestr(newtime,'yyyymmdd-HHMMSS');
                    newtime2=newtime2(end-5:end-4);
                %Compute the new juliand day
                if day(newtime,'dayofyear')<10
                    juliandnew=['00' num2str(day(newtime,'dayofyear'))];
                elseif day(newtime,'dayofyear')>=10 & day(newtime,'dayofyear')<100
                    juliandnew=['0' num2str(day(newtime,'dayofyear'))]; 
                elseif day(newtime,'dayofyear')>=100 
                    juliandnew=num2str(day(newtime,'dayofyear')); 
                end
                %Create the new name of the signal file
                newfile=strcat(mainpath,'/',  juliandnew,'/',sta_ref,'.BH',comp,'.',num2str(year(newtime)),'.',juliandnew,'.',newtime2,'0000.sac');
                %Control the directory of the new signal file 
                p3=dir(newfile);
                if length(p3)~=0
                    %Load the trace
                    load(newfile); 
                    z3=mytrace.data;nt=mytrace.sampleCount;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate;
                    %Rescale and remove the best straight-fit line
                    z3=z3-mean(z3);
                    z3=detrend(z3);
                    %Extraction of the waveform
                    Z1=z(TP_sampl-w1:end);
                    Z2=z3(1:TP_sampl+w1-length(z));
                    Z=[Z1;Z2];
                else
                    %Warning message
                    set(text1, 'String', 'Warning. One or more signals are no available.');
                    drawnow
                end
                elseif TP_sampl-w1<0 & TP_sampl+w2>length(z)
                    %Warning message
                    set(text1, 'String', 'Warning. One or more signals are no available.');
                    drawnow 
                    continue 
                end
                %Control the length of the trace
                if length(Z)==round(fs)*2*wlong+1
                    %Create the Output filename and the directory 
                    timeref=datestr(Timesel,'yyyymmdd-HHMMSS');
                    timeref=timeref(end-5:end);
                    output=strcat(pathout,juliand,'/',sta_ref,'.BH',comp,'.',Name(end-18:end-10),timeref,outputformat);
                    mytrace.data=Z;mytrace.sampleCount=length(Z);mytrace.startTime=datenum(Timesel-seconds(wlong));mytrace.endTime=datenum(Timesel+seconds(wlong));
                    if ~exist(strcat(pathout,juliand,'/'), 'dir');mkdir(strcat(pathout,juliand,'/'));end
                    %Save the extracted waveform for every station 
                    if outputformat=='.mat'    
                        save (output,'mytrace');
                    elseif outputformat=='.txt' 
                        Time=datenum([Timesel-seconds(wlong):seconds(dt):Timesel+seconds(wlong)]);
                        data=table(Time',Z);data.Properties.VariableNames={'Time', 'Amplitude'};
                        data.Properties.VariableUnits = {'datenum','counts'};
                        writetable(data,output)
                    end
                end
                %Command message
                set(text1, 'String', sprintf(['Signal-' num2str(ii) ' of ' num2str(bb) outputformat '-Observed data are saved.']));
                drawnow
            else
                %Error message
                set(text1, 'String', sprintf(['Error. Station names are no found.']));
                drawnow
                return
            end
        end
    end
catch
    %Error assessment
    if isempty(Timesel);set(text1, 'String','Error. Invalid trigger value.');drawnow;return;end
    if isempty(Name);set(text1, 'String','Error. Trigger matrix is empty.');drawnow;return;end 
    if isempty(Z);set(text1, 'String','Error. No correct signals data.');drawnow;return;end  
end
end