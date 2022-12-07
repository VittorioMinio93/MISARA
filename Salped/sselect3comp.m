%%%%%%%%Select the three components of  the waveforms detected%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sta_ref: the name of the reference station
% comp, component: the name of the components
% Timesel: trigger time (datenum format)
% Name: the path of the file of reference 
% mainpath: the Input folder
% fs: sampling rate
% dt: sampling interval
% nt: number of samples
% startTime, endTime: temporal limits of the traces 
% wshort, wlong: the pre-trigger and post-trigger window (s)
% pathout: the Output folder
% outputformat:  the format (.mat/.txt) of the output files 
% text1:  Window message 
% z,z3: amplitude vector of the seismic trace
% Z,ZZ: one component waveform extracted
% YY: three components waveforms extracted
% Zz,Nn,Ee: the three components
% mytrace: seismic trace structure
% TP_sampl: trigger time (samples)
% juliand, juliandnew: julian day

function [Zz,Nn,Ee]=sselect3comp(Timesel,Name,mainpath,fs,wlong,wshort,sta_ref,text1)
try 
    ZZ=[];%Initialize the waveform extracted vector 
    %Set the component names
    components={'Z','N','E'};
    YY=[];%Initialize the three components waveforms extracted 
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
        %Create the name of the signal file 
        file_sac=strcat(mainpath,'/',juliand,'/',sta_ref,'.BH',comp,Name(end-19:end));
        list=dir(file_sac);
        ZZ=[];%Initialize the waveform extracted vector 
        for kk=1:1
            %Control the directory of the signal file
            filesac=strcat(mainpath,'/',juliand,'/',list(kk).name);
            p=dir(filesac);
            if ~isempty(p)
                %Load the trace
                load(filesac); 
                z=mytrace.data;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate; nyq=0.5/dt; nt=mytrace.sampleCount;
                startTime=mytrace.startTime;endTime=mytrace.endTime;
                %Rescale and remove the best straight-fit line
                z=z-mean(z);
                z=detrend(z);
                w1=round(wshort*fs);
                w2=round(wlong*fs);
                Z=[];%Initialize the waveform extracted vector 
                %Extracting the waveforms on the basis of the trigger time  and signal files available
                if TP_sampl-w1>0 & TP_sampl+w2<=length(z)
                    %Extraction of the waveform
                    Z=z(TP_sampl-w1:TP_sampl+w2);
                elseif TP_sampl-w1<0 & TP_sampl+w2<=length(z)
                    %Search the new time of reference
                    newtime=Timesel-1/24;
                    newtime2=datestr(newtime,'yyyymmdd-HHMMSS');
                    newtime2=newtime2(end-5:end-4);
                    %Compute the new juliand day
                    if day(newtime,'dayofyear')<10
                        juliandnew=['00' num2str(day(newtime,'dayofyear'))]
                    elseif day(newtime,'dayofyear')>=10 & day(newtime,'dayofyear')<100
                        juliandnew=['0' num2str(day(newtime,'dayofyear'))]; 
                    elseif day(newtime,'dayofyear')>=100 
                        juliandnew=num2str(day(newtime,'dayofyear')); 
                    end
                    %Create the new name of the signal file
                    newfile=strcat(path,'/',juliandnew,'/',sta_ref,'.BH',comp,'.',num2str(year(newtime)),'.',juliandnew,'.',newtime2,'0000.sac');
                    %Control the directory of the new signal file 
                    p3=dir(newfile);
                    if ~isempty(p3)
                        %Load the trace
                        load(newfile); 
                        z3=mytrace.data;nt=mytrace.sampleCount;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate;
                        %Rescale and remove the best straight-fit line
                        z3=z3-mean(z3);
                        z3=detrend(z3);
                        %Extraction of the waveform
                        Z1=z(1:TP_sampl+w2);
                        Z2=z3(end-(w1-TP_sampl):end);
                        Z=[Z2;Z1];
                    else
                        %Warning message
                        set(text1, 'String', 'Warning. One or more signals are no available.');
                        drawnow
                    end
                elseif TP_sampl-w1>0 & TP_sampl+w2>length(z)
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
                    newfile=strcat(path,'/',juliandnew,'/',sta_ref,'.BH',comp,'.',num2str(year(newtime)),'.',juliandnew,'.',newtime2,'0000.sac');
                    %Control the directory of the new signal file
                    p3=dir(newfile);
                    if ~isempty(p3)
                        %Load the trace
                        load(newfile); 
                        z3=mytrace.data;nt=mytrace.sampleCount;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate;
                        %Rescale and remove the best straight-fit line
                        z3=z3-mean(z3);
                        z3=detrend(z3);
                        %Extraction of the waveform
                        Z1=z(TP_sampl-w1:end);
                        Z2=z3(1:TP_sampl+w2-length(z));
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
                if length(Z)==round(fs)*(wshort+wlong)+1
                    ZZ=[ZZ Z];
                end
            else
                %Error 
                return
            end            
        end
        %Concatenate through the columns
        YY=[YY ZZ];
    end
    %Control the size of the three components matrix 
    if size(YY,2)==3
        Zz=YY(:,1); Nn=YY(:,2);Ee=YY(:,3);
    else
        %Error message
        set(text1, 'String','Error. One or more components are no available.')
        drawnow
        return               
    end
catch
    %Error assessment
    if isempty(Timesel);set(text1, 'String','Error. Invalid trigger value.');drawnow;return;end
    if isempty(Name);set(text1, 'String','Error. Trigger matrix is empty.');drawnow;return;end 
    if isempty(list);set(text1, 'String','Error. No files are found.');drawnow;return;end
    if isempty(ZZ);set(text1, 'String','Error. No correct signals data.');drawnow;return;end    
end
end