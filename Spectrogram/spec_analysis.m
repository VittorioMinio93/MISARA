%%%%%%%%%%%%% Calculation of the Spectrogram vs. Rms %%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params: user parameters
% avfig: Spectrogram figure
% panel_color: colors of the main panel 
% fig_color: colors of the figure
% pan1: control right side panel 
% text1: Window message
% axes1,axes3: axes of the current figure
% f1_lim, f2_lim: limit of y-axis on spectrogram axes
% t1_lim, t2_lim: limit of x-axis on spectrogram and rms axes
% slidval: value of dynamic plots 
% logval: set the scale of the spectral amplitude
% specmeanval: Active/Deactive the averaging of the spectrums
% h: UI control parameters
%%%%%%%%%%%h.variables/params.variables/variables%%%%%%%%%%%%%%%%%%%%%%%%%%
% rad: radio button within a button group 
% rad1, rad2, rad3, rad4: button group 
%-------------------------------- General settings
%--- General info
% sta_ref: the name of the reference station 
% component, comp: the name of the component/channel
% comp_chan: the station system used (component/channel)
%--- Data file
% data: the Input folder
% save: the Output folder
% outputformat: the format (.mat/.txt) of the output files
% search_sac: term of research file 
% pathname, fname: Input folder and list of filenames
% ff: list of files with their path
% mytrace: seismic trace structure
% dt: sampling interval
% fs: sampling rate
% nyq: nyquist number
% z: amplitude vector of the seismic trace
% nt,NT: number of samples
% startTime, endTime: temporal limits of the traces  
%--- Spectral setting
% w_spec: the analysis window for the Spectrogram (s)
% np: the dimension of the spectrums
% f_high_pass: the frequency for the high-pass filter (Hz)
% spec_hour, numspec: the number of spectrums for each trace
% Time_spec, TIME_SPEC: time vectors of the spectrogram 
% S, SSTOT: spectrum vectors
%--- RMS setting
% w_rms: the analysis window for the RMS analysis (s)
% f1: the minimum frequency for the RMS analysis (Hz)
% f2: the maximum frequency for the RMS analysis (Hz)
% step_rms: the frequency step for the RMS analysis (Hz)
% factor: the temporal factor for the moving average of the RMS (hours) 
% Time_rms, TIME_RMS, TIME_AMP: time vectors of the rms
% amp, AMP,  AMPP: amplitude vectors of the rms
% F:  range of frequency analysis (Hz) 

function spec_analysis()
%% Declaration global parameters
clc
global params
comp=params.component;
component=[];
ff= cell(1);
mytrace=[];
numspec=[];
TIME_AMP=[];
TIME_SPEC=[];
AMPP=[];
SSTOT=[];
slidval=1;
dt=[];
fs=[];
nyq=[];
F=[];
f1=[];f2=[];step_rms=[];factor=[];w_rms=[];
f1_lim=[];f2_lim=[];
t1_lim=[];t2_lim=[];
logval=0;
specmeanval=1;

%% Spectrogram figure
fig_color= [.94 .94 .94];
panel_color= [.97 .97 .97];

avfig= figure('Name', 'Spectrogram', 'position',get(0,'screensize'), ...
    'NumberTitle','off','units','normalized', 'color', fig_color,'MenuBar','none');%%Create the spectrogram figure 

%Window message to user
text1= uicontrol('Style','edit','FontName','Arial','Fontsize',10,'Fontweight','bold',...
                'String','Window Message','Max',5,'Min',0,'HorizontalAlignment','left',...
                'BackgroundColor', 'w','Enable','inactive',...
                'units','normalized','Position',[.05 .02 .72 .05]);


%% Axeses
%Spectrogram axes
axes1= axes('position',[.05 .56 .72 .40]);
hold(axes1,'on')
colorbar(axes1,'northoutside')
colormap(axes1,'parula')
ylabel(axes1,'Frequency(Hz)');
ylim(axes1,[0 15]);
% caxis(axes1,[0 1])

%Rms axes
axes3= axes('position',[.05 .14 .72 .36]);
hold(axes3,'on')
grid(axes3,'on')
grid(axes3,'minor')
ylabel(axes3,'Rms (a.u.)');
ylim(axes3,[0 1])

%% Control right side panel
pan1= uipanel(avfig,'visible','on','Position',[.81 .03 .185 .94],...
    'BackgroundColor', panel_color);

uipanel(avfig,'visible','on','Position',[.81 .45 .185 .01],...
    'BackgroundColor', fig_color);

%% Text boxes
%-------------------------------- Axis setting
uicontrol(pan1,'Style','text', 'String','Axis setting',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.03 .73 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','F1 (Hz)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .70 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','F2 (Hz)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .65 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','t1',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .58 .2 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','t2',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .53 .2 .04],...
    'BackgroundColor',panel_color);

%-------------------------------- RMS setting
uicontrol(pan1,'Style','text', 'String','Rms setting',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.55 .73 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Win (s)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .70 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','F1 (Hz)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .65 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','F2 (Hz)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .60 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Step (Hz)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .53 .2 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Factor (h)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .48 .4 .04],...
    'BackgroundColor',panel_color);

%% Editable text 
%-------------------------------- Axis setting
h.f1_lim= uicontrol(pan1,'Style','edit', 'String','0',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.2 .69 .3 .04],...
    'BackgroundColor','w','callback',@setlimits1);
h.f2_lim= uicontrol(pan1,'Style','edit', 'String','15',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.2 .64 .3 .04],...
    'BackgroundColor','w','callback',@setlimits1);
h.t1_lim= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.2 .59 .3 .04],...
    'BackgroundColor','w','callback',@setlimits2);
h.t2_lim= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.2 .54 .3 .04],...
    'BackgroundColor','w','callback',@setlimits2);

%-------------------------------- RMS setting
h.w_rms= uicontrol(pan1,'Style','edit', 'String',params.w_rms,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .69 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.f1= uicontrol(pan1,'Style','edit', 'String',params.f1,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .64 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.f2= uicontrol(pan1,'Style','edit', 'String',params.f2,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .59 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.step_rms= uicontrol(pan1,'Style','edit', 'String',params.step_rms,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .54 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.factor= uicontrol(pan1,'Style','edit', 'String',params.factor,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .49 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);

%% Add Radio buttons
% Editable text 
switch params.comp_chan
    case'Comp'
     uicontrol(pan1,'Style','text', 'String','Component',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.05 .85 .4 .05],...
    'BackgroundColor',panel_color);
    case'Chan'
    uicontrol(pan1,'Style','text', 'String','Channel',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.05 .85 .4 .05],...
    'BackgroundColor',panel_color);
end

% Component selection radio buttons
h.rad= uibuttongroup(pan1,'units','normalized','BackgroundColor',panel_color,...
    'bordertype','none','Position',[.07 .82 .9 .05]);

set(h.rad,'SelectionChangeFcn',@radcbk);

h.rad1 = uicontrol( h.rad, 'Style','Radio','String','Z',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.0 0 .15 1],'HandleVisibility','off');
h.rad2 = uicontrol( h.rad, 'Style','Radio','String','N',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.15 0 .15 1],'HandleVisibility','off');
h.rad3 = uicontrol( h.rad, 'Style','Radio','String','E',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.3 0 .15 1],'HandleVisibility','off');
h.rad4 = uicontrol( h.rad, 'Style','Radio','String','F',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.45 0 .15 1],'HandleVisibility','off');

%Set the component/channel of reference 
if strcmp(comp,'Z')
    set (h.rad1,'value',1);
    set (h.rad2,'value',0);
    set (h.rad3,'value',0);
    set (h.rad4,'value',0);

elseif strcmp(comp,'N')
    set (h.rad1,'value',0);
    set (h.rad2,'value',1);
    set (h.rad3,'value',0);
    set (h.rad4,'value',0);

elseif strcmp(comp,'E')
    set (h.rad1,'value',0);
    set (h.rad2,'value',0);
    set (h.rad3,'value',1);
    set (h.rad4,'value',0);
elseif strcmp(comp,'F')
    set (h.rad1,'value',0);
    set (h.rad2,'value',0);
    set (h.rad3,'value',0);
    set (h.rad4,'value',1);
end

%% Buttons
% Open files button
uicontrol('style','pushbutton', 'string','Open files','units','normalized',...
    'position',[.815 .9 .08 .0494], 'callback',@openfnc);


% Open folder button
uicontrol('style','pushbutton', 'string','Open folder','units','normalized',...
    'position',[.91 .9 .08 .0494], 'callback',@openfolder);

% Editable text 
uicontrol(pan1,'Style','text', 'String','Avg Mode',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.66 .85 .4 .05],...
    'BackgroundColor',panel_color);
% Averaging spectrums checkbox
uicontrol('Style','checkbox','BackgroundColor', panel_color,...
                'String','on/off','units','normalized','Tag','logcb',...
                'Value',1,'Position',[.94 .805 .05 .03],'callback',@specmean);

% Logarithmic scale checkbox
uicontrol('Style','checkbox','BackgroundColor', panel_color,...
                'String','Log Axis','units','normalized','Tag','logcb',...
                'Value',0,'Position',[.815 .49 .05 .03],'callback',@logscale);
% Logarithmic spectrum checkbox
uicontrol('Style','checkbox','BackgroundColor', panel_color,...
                'String','Log Spec','units','normalized','Tag','logcb',...
                'Value',0,'Position',[.815 .46 .05 .03],'callback',@logspec);
% Editable texts
uicontrol(pan1,'Style','text', 'String','Spectrogram',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.03 .33 .4 .07],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Rms',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.55 .33 .4 .07],...
    'BackgroundColor',panel_color);
% Calculation and save button
uicontrol('style','pushbutton', 'string','Calculation','units','normalized',...
    'position',[.815 .3 .08 .0494], 'callback',@calc); 
uicontrol('style','pushbutton', 'string','Calculation','units','normalized',...
    'position',[.91 .3 .08 .0494], 'callback',@calc2);
uicontrol('style','pushbutton', 'string','Save Plot','units','normalized',...
    'position',[.815 .2 .08 .0494], 'callback',@save1); 
uicontrol('style','pushbutton', 'string','Save Plot','units','normalized',...
    'position',[.91 .2 .08 .0494], 'callback',@save2);
uicontrol('style','pushbutton', 'string','Save Data','units','normalized',...
    'position',[.815 .1 .08 .0494], 'callback',@savespec);
uicontrol('style','pushbutton', 'string','Save Data','units','normalized',...
    'position',[.91 .1 .08 .0494], 'callback',@saverms);

% Next and pre button
uicontrol('style','pushbutton', 'string','Next','units','normalized',...
    'position',[.776 .4 .03 .04], 'callback',@rms2);
uicontrol('style','pushbutton', 'string','Pre','units','normalized',...
    'position',[.776 .46 .03 .04], 'callback',@rms1);

%% Other functions (nested)
%%%%Return the list of files manually selected%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function openfnc(~,~)
    TIME_AMP=[];TIME_SPEC=[];AMPP=[];SSTOT=[];%Initialize the spectrogram and rms vectors
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);cla(axes3);%Clear the axes
    tt= 0;
    ff=cell(1);%Initialize the list of files
    search_sac=strcat(params.sta_ref,'.BH',comp,'*.mat');
    component=comp;
    [fname, pathname] = uigetfile({search_sac;'*.*'},'File Selector',...
        params.data,'MultiSelect','on');%Selection of the files    
    
    [~,nnn]= size(fname);
    if nnn==1 && fname==0
        %Error message
        set(text1, 'String', 'Error. No files selection made.');
        drawnow        
    else        
        fname= cellstr(fname);
        [~,np]= size(fname);%Number of files selected
        if np<2
            ff=cell(1);%Initialize the list of files
            tt= 1;
            %Error message
            textLabel = sprintf('Error. Not sufficient number of files.');
            set(text1, 'String', textLabel);
            drawnow
        else
            %Loop through the number of files selected
            for ii=1:np
                ff{ii}= [pathname,fname{ii}];%Return the list of files with their paths
            end
            %Command message
            textLabel = sprintf('Selected files are opened.');
            set(text1, 'String', textLabel);
            drawnow
        end                
    end
end

%==========================================================================
%%%%Return the list of files selected in a specific folder%%%%%%%%%%%%%%%%%
function openfolder(~,~)
    TIME_AMP=[];TIME_SPEC=[];AMPP=[];SSTOT=[];%Initialize the spectrogram and rms vectors
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);cla(axes3);%Clear the axes
    pathname= uigetdir(params.data);
    search_sac=strcat(params.sta_ref,'.BH',comp,'*.mat');
    component=comp;
    pathname=[pathname '/'];
    if pathname~=0
    % Get all files in the directory and subdirectories
    [ff,~]= getAllFiles(pathname,search_sac);
    ff=ff';
    
    if ~isempty(ff) 
        %Command message
        textLabel = sprintf('All files in the folder are opened.');
        set(text1, 'String', textLabel); 
        drawnow
    else
        ff=cell(1);%Initialize the list of files
        %Error message
        textLabel = sprintf('Error. No files are read.');
        set(text1, 'String', textLabel);
        drawnow
    end
    
    else
         ff=cell(1);%Initialize the list of files
        %Error message
        set(text1, 'String', 'Error. No folder selection made.');
        drawnow
    end
end

%==========================================================================
%%%%Return the name of the component/channel%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radcbk(~,eventdata)
    
    comp = get(eventdata.NewValue,'String');   
    
end

%==========================================================================
%%%%Calculation and plot of the spectrogram%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calc(~,~) 
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1)%Clear the axes
    try
      %Validation of the analysis parameters
      if params.w_spec<=0; set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
      if params.np<=0; set(text1, 'String','Error. Invalid Spectrogram parameters.');drawnow;return;end
      if params.f_high_pass<=0;set(text1, 'String','Error. Invalid High pass filter frequency.');drawnow;return;end
      if params.spec_hour<=0; set(text1, 'String','Error. Invalid Spectrogram parameters.');drawnow;return;end
      %Command message
      textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
      set(text1, 'String', textLabel); 
      drawnow
      %Reading .mat files
      [~, len]=size (ff);%Number of files selected
      SSTOT=cell(1,len);%Initialize the spectrum matrix
      NT=0;
      TIME_SPEC=cell(len,1);%Initialize the time vector
      %Loop through the number of files selected
      for i=1:len
          load(ff{i});%load the i-th file
          z=mytrace.data;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate; nyq=0.5/dt; nt=mytrace.sampleCount;
          startTime=mytrace.startTime;endTime=mytrace.endTime;
          %Rescale and remove the best straight-fit line
          z=z-mean(z);
          z=detrend(z);  
          %Calculation of the spectrum and time vectors
          [S,numspec]=SPEC(z,dt,fs,nyq,nt);
          [Time_spec,~]=TIME(startTime,endTime,nt);
          %Concatenate the results
          SSTOT(i)={S};
          NT=NT+nt;
          TIME_SPEC(i)={datenum(unique(Time_spec))};
      end
      %Plot the results in the spectrogram axes
      ploting1(TIME_SPEC,SSTOT,nyq,params.np,axes1)
      %Set the limits of x-axis on the spectrogram and rms axes
      time=datetime(cell2mat(TIME_SPEC),'ConvertFrom','datenum');%Convert in datetime format
      if ~isempty(time) & ~isempty(cell2mat(TIME_AMP))
          xlim(axes1,[min(time) max(time)])
          xlim(axes3,[min(time) max(time)])
      elseif ~isempty(time) & isempty(cell2mat(TIME_AMP)) 
          xlim(axes1,[min(time) max(time)])
      end
      %Command message
      set(text1, 'String','Calculation finished.');
      drawnow
    catch
        %Error assessment
        if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
        if isempty(z);set(text1, 'String','Error. No correct signals data.');drawnow;return;end
        if length(z)<=round(params.w_spec*fs);set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
    end    
end
%==========================================================================
%%%%Calculation and plot of the rms%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calc2(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes3)%Clear the axes
    try
        %Acquisition and validation of the analysis parameters
        w_rms=str2num(get(h.w_rms,'string'));if w_rms<=0;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
        factor=str2num(get(h.factor,'string'));if factor<=0;set(text1, 'String','Error. Invalid temporal  factor mean.');drawnow;return;end   
        f1=str2num(get(h.f1,'string'));if f1<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        f2=str2num(get(h.f2,'string'));if f2<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        step_rms=str2num(get(h.step_rms,'string'));if step_rms<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        F=f1:step_rms:f2;
        if isempty(F)| length(F)<2;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        %Command message
        textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
        set(text1, 'String', textLabel); 
        drawnow
        %Reading .mat files
        [~, len]=size (ff);%Number of files selected
        TIME_AMP=cell(len,length(F)-1);%Initialize the time matrix
        AMPP=cell(len,length(F)-1);%Initialize the rms matrix
        %Loop through the frequency analysis
        for k=1:length(F)-1
            TIME_RMS=cell(len,1);%Initialize the time vector
            AMP=cell(len,1);%Initialize the rms vector
            %Loop through the number of files selected
            for i=1:len
                load(ff{i});%load the i-th file 
                z=mytrace.data;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate; nyq=0.5/dt; nt=mytrace.sampleCount;
                startTime=mytrace.startTime;endTime=mytrace.endTime;
                %Rescale and remove the best straight-fit line
                z=z-mean(z);
                z=detrend(z);
                %Calculation of the rms and time vectors
                amp=RMS(z,fs,dt,F(k),F(k+1));
                [~,Time_rms]=TIME(startTime,endTime,nt);
                %Concatenate the results
                TIME_RMS(i)={datenum(Time_rms)};
                AMP(i)={amp};
            end
            %Concatenate the results through the columns
            TIME_AMP(:,k)= TIME_RMS;
            AMPP(:,k)= AMP;
        end
        %Set the block of results
        slidval=length(F)-1;
        %Plot the results in the rms axes
        ploting2(slidval,axes3)
        %Set the limits of x-axis on rms axes
        if ~isempty(cell2mat(TIME_SPEC))
            time=datetime(cell2mat(TIME_SPEC),'ConvertFrom','datenum');%Convert in datetime format
            xlim(axes3,[min(time) max(time)])
        end
       %Command message
       set(text1, 'String','Calculation finished.');
       drawnow
    catch
        %Error assessment
        if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
        if isempty(z);set(text1, 'String','Error. No correct signals data.');drawnow;return;end
        if length(z)<=round(w_rms*fs);set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
    end
end

%==========================================================================
%%%%Calculation of the time vector%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [time_spec,time_rms]=TIME(startTime,endTime,nt)
    time_spec=[];%Initialize the time vector of the spectrogram
    time_rms=[];%Initialize the time vector of the rms
    %Convert the time limits of the trace in datetime format
    Time1=datetime(startTime,'ConvertFrom','datenum');
    Time2=datetime(endTime,'ConvertFrom','datenum');
    
    if specmeanval==1
    %Compute the time vector of the spectrums
    Time_spec_1=Time1+numspec*(seconds(params.w_spec+1/fs));
    Time_spec_2=Time1+params.spec_hour*numspec*(seconds(params.w_spec+1/fs));
    time_spec=[Time_spec_1:numspec*(seconds(params.w_spec+1/fs)):Time_spec_2]';
    else
    Time_spec_1=Time1+seconds(params.w_spec);
    Time_spec_2=Time2;
    Time_spec=Time_spec_1:seconds(params.w_spec):Time_spec_2;
    time_spec=[time_spec;Time_spec'];
        
    end
    %Compute the time vector of the rms
    Time_rms_1=Time1+seconds(w_rms);
    Time_rms_2=Time2;
    Time_rms=Time_rms_1:seconds(w_rms):Time_rms_2;
    time_rms=[time_rms;Time_rms'];   
end
%================================
%%%%Calculation of the spectrums %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [STOT,numspec]=SPEC(z,dt,fs,nyq,nt)
    try
        STOT=[];%Initialize the spectrogram vector
        %Filter the seismic traces with a high-pass
        ff1=params.f_high_pass;
        fn=fs/2;
        W=ff1/fn;
        [x y]=butter(2,W,'high');
        z=filter(x,y,z);
        SS=[];%Initialize the spectrogram vector
        jj=0;i1=0;i2=0;
        %Move the analysis window through the seismic traces
        wl=round(params.w_spec/dt);np=params.np;
         while i2 < nt-wl || i2== nt-wl
             jj=jj+1;
             i1=(1+(jj-1)*wl); i2=i1+wl-1;%Limits of the analysis window
             hwindow=hann(length(z(i1:i2)));%Returns a symmetric Hann window
             %Apply the Fast Fourier Transform algorithm and compute the spectrums
             spec=fft(z(i1:i2).*hwindow,np);
             spec=spec.*conj(spec);
             spec=spec(1:np/2);
             %Concatenate the results
             SS=[SS,spec];
         end
         
         %Averaging Mode on/off
         if specmeanval==1
         %Meaning of the spectrums based on the numspec parameter
         [nrow nspec]=size(SS);%Dimensions of the spectrums matrix
         nseconds=round(nt/fs);
         numspec=(nseconds/params.w_spec)/params.spec_hour;%Number of the spectrums per seismic trace 
         %Loop through the number of spectrums computed
         for kk=1:numspec:nspec
             ii1=(1+(kk-1)); ii2=ii1+numspec-1;
             SPE=mean(SS(:,ii1:ii2)');
             %Concatenate the results
             STOT=[STOT,SPE'];
         end
         else
         %Loop through the number of spectrums computed
         [nrow nspec]=size(SS);%Dimensions of the spectrums matrix
         numspec=0; 
         for kk=1:nspec
             % Calculate spectrogram in dB
             SPE= log10(abs(SS(:,kk)));
             %Concatenate the results
             STOT=[STOT,SPE];
         end
         end
    catch
        %Error assessment
        if params.f_high_pass>fn;set(text1, 'String','Error. Invalid High pass filter frequency.');drawnow;return;end
    end
end

%==========================================================================
%%%%Calculation of the rms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function amp=RMS(z,fs,dt,f1,f2)
    try
        %Filter the seismic traces with a band-pass
        fs=(1/dt);fn=fs/2;
        W(1,1)=f1/fn; W(1,2)=f2/fn;
        [a b]=butter(2,W);
        z=filter(a,b,z);
        %Move the analysis window through the seismic traces
        N=length(z);wl=round(w_rms/dt);
        i1=0; i2=0;k=0;
        while i2 < N-wl || i2== N-wl
            k=k+1;
            i1=(1+(k-1)*wl); i2=i1+wl-1;%Limits of the analysis window
            y1=z(i1:i2).*tukeywin(wl);%tapered cosine window
            %Apply the rms algorithm
            amp(k)=sqrt(sum((y1.*conj(y1))/wl));
        end
        amp=amp';
    catch
        %Error assessment
        if f1>=fn | f2>fn | step_rms>=fn;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    end
end

%==========================================================================
%%%%Plot of the spectrogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting1(time,spec,nyq,np,ax)
    try
        Freq=[];%Initialize the frequency vector
        Time=[];%Initialize the time vector
        cla(ax)%Clear the spectrogram axes
        time=datetime(cell2mat(time),'ConvertFrom','datenum');%Convert in datetime format
        spec=cell2mat(spec);
        [Time Freq]=meshgrid(time,(nyq/(np/2)):(nyq/(np/2)):nyq);
        spec=spec(1:size(Freq,1),:);%Select the spectrums based on the nyquist number
        
        if specmeanval==1
        %Normalize the spectrums
        for ii=1:1:size(spec,2)
            spec(:,ii)=spec(:,ii)./max(spec(:,ii));
        end
        c1=0;c2=1;
        else
        c1=prctile(spec(:),40);
        c2=prctile(spec(:),99.9);
        end
        %Logarithmic spectrums
        if logval==1
            spec=log10(spec);
        end
        %Plot the spectrogram
        surf(ax,Time,Freq,spec)
        view (2)
        shading(ax,'interp')
        axis(ax,'tight')
        colorbar(ax,'northoutside')
        colormap(ax,'parula')
        ylabel(ax,'Frequency(Hz)');
        ylim(ax,[0 15]);
        caxis(ax,[c1,c2])
    catch
        %Error assessment
        if size(Time)~=size(Freq) & size(Time)~=size(spec) & size(Freq)~=size(spec)
            set(text1, 'String','Error. Incoherent dimensions of the Spectrogram.');drawnow;return;
        end       
    end
end

%==========================================================================
%%%%Plot of the rms at different frequency analysis%%%%%%%%%%%%%%%%%%%%%%%%
function ploting2(ii,ax)
    try
        cla(ax)%Clear the rms axes
        amp=cell2mat(AMPP);
        time=datetime(cell2mat(TIME_AMP),'ConvertFrom','datenum');%Convert in datetime format
        amp=amp(:,ii);%Select the ii-th block of the rms
        time=time(:,ii);%Select the ii-th block of the time
        rms_mean=movmean(amp,hours(factor),'SamplePoints',time);%Moving average of the rms according the factor value
        %Plot the normalized amplitude of the rms
        plot(ax,time,rms_mean./max(rms_mean))
        legend(ax,[num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
        ylabel(ax,'Rms (a.u.)');
        ylim(ax,[0 1])
    catch
        %Error assessment
        if size(amp,1)~=size(time,1);set(text1, 'String','Error. Incoherent dimensions of the Rms array.');drawnow;return;end
    end
end

%==========================================================================
%%%%Set the scale of the y-axis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logscale(h, ~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    if (get(h,'Value') == get(h,'Max'))
        %Command message
        set(text1, 'String', 'Logarithmic axis chosen.');
        %Set the semi-logarithmic scale
        set(axes3,'yScale','log'); 
        set(axes1,'yScale','log');
    else
        %Command message
        set(text1, 'String', 'Linear axis chosen.');
        %Set the linear scale
        set(axes1,'yScale','linear');
        set(axes3,'yScale','linear');
    end
end

%==========================================================================
%%%%Set the scale of the spectral amplitude%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logspec(h, ~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Error control
   if isempty(cell2mat(SSTOT));set(text1, 'String','Error. Spectrogram matrix is empty.');drawnow;return;end
    if (get(h,'Value') == get(h,'Max'))
        %Command message
        set(text1, 'String', 'Logarithmic spectrum values.');
        %Set the logarithmic values of the spectrums
        logval=1;
        ploting1(TIME_SPEC,SSTOT,nyq,params.np,axes1)
    else
        %Command message
        set(text1, 'String', 'Linear spectrum values.');
        %Set the logarithmic values of the spectrums
        logval=0;
        ploting1(TIME_SPEC,SSTOT,nyq,params.np,axes1)
    end
end

%==========================================================================
%%%%Active/Deactive the averaging of spectrums%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function specmean(h, ~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    if (get(h,'Value') == get(h,'Max'))
        %Command message
        set(text1, 'String', 'Averaging mode is Activated.');
        %Set the average computing of the spectrums
        specmeanval=1;
    else
        %Command message
        set(text1, 'String', 'Averaging mode is Deactivated.');
        %Set the computing of the spectrums
        specmeanval=0;
    end
end
%==========================================================================
%%%%Refresh the rms axes through the pre button%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rms1(~,~)
   %Error control
   if isempty(cell2mat(AMPP));set(text1, 'String','Error. Rms matrix is empty.');drawnow;return;end
   %Decrease the slidval value
   slidval=slidval-1;if slidval<=0;slidval=1;end
   %Plot the new block of rms 
   ploting2(slidval,axes3);        
end

%==========================================================================
%%%%Refresh the rms axes through the next button%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rms2(~,~)
    %Error control
    if isempty(cell2mat(AMPP));set(text1, 'String','Error. Rms matrix is empty.');drawnow;return;end
    %Increase the slidval value
    slidval=slidval+1;if slidval>length(F)-1;slidval=length(F)-1;end
    %Plot the new block of rms 
    ploting2(slidval,axes3);        
end

%==========================================================================
%%%%Set the limits of the y-axis on the spectrogram axes%%%%%%%%%%%%%%%%%%%
function setlimits1(~,~)
    try
        %Manual settings
        %Acquisition frequency limits
        f1_lim=str2num(get(h.f1_lim,'string'));
        f2_lim=str2num(get(h.f2_lim,'string')); 
        ylim(axes1,[f1_lim f2_lim])
    catch
        %Default settings
        ylim(axes1,[0 15])
        %Reset frequency limits
        set(h.f1_lim,'string','0');
        set(h.f2_lim,'string','15');      
    end
end

%==========================================================================
%%%%Set the limits of the x-axis on the spectrogram and rms axes%%%%%%%%%%%
function setlimits2(~,~)
    try
        %Manual settings
        %Acquisition time limits
        t1_lim=datetime(get(h.t1_lim,'string'),'InputFormat','yyMMdd-HHmmSS');%Convert in datetime format
        t2_lim=datetime(get(h.t2_lim,'string'),'InputFormat','yyMMdd-HHmmSS'); %Convert in datetime format
        xlim(axes1,[t1_lim t2_lim])
        xlim(axes3,[t1_lim t2_lim])
    catch
        %Default settings
        if ~isempty(cell2mat(TIME_SPEC)) 
            time=datetime(cell2mat(TIME_SPEC),'ConvertFrom','datenum');%Convert in datetime format
            xlim(axes1,[min(time) max(time)])
            xlim(axes3,[min(time) max(time)])
        end
    end
end
%==========================================================================
%%%%Control and reset automatically analysis parameters%%%%%%%%%%%%%%%%%%%%
function deflimits(~,~)
    %Acquisition analysis parameters
    w_rms=str2num(get(h.w_rms,'string'));
    factor=str2num(get(h.factor,'string'));   
    f1=str2num(get(h.f1,'string'));
    f2=str2num(get(h.f2,'string'));
    step_rms=str2num(get(h.step_rms,'string'));
    
    %Control and eventually reset parameters
    if (size(w_rms,1)==0);set(h.w_rms,'string',params.w_rms);end
    if (size(factor,1)==0);set(h.factor,'string',params.factor);end
    if (size(f1,1)==0);set(h.f1,'string',params.f1);end
    if (size(f2,1)==0);set(h.f2,'string',params.f2);end
    if (size(step_rms,1)==0);set(h.step_rms,'string',params.step_rms);end
end

%==========================================================================
%%%%Save the spectrogram plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save1(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    try
        Freq=[];%Initialize the frequency vector
        Time=[];%Initialize the time vector
        time=datetime(cell2mat(TIME_SPEC),'ConvertFrom','datenum');%Convert in datetime format
        spec=cell2mat(SSTOT);
        [Time Freq]=meshgrid(time,(nyq/(params.np/2)):(nyq/(params.np/2)):nyq);
        spec= spec(1:size(Freq,1),:);%Select the spectrums based on the nyquist number
        if specmeanval==1        
        %Normalize the spectrums
        for ii=1:1:size(spec,2)
            spec(:,ii)= spec(:,ii)./max( spec(:,ii));
        end
        c1=0;c2=1;
        else
        c1=prctile(spec(:),40);
        c2=prctile(spec(:),99.9);
        end
        %Logarithmic spectrums
        if logval==1
            spec=log10(spec);
        end
        %Plot the spectrogram
        new_figure=figure('Name','Spectrogram plot','NumberTitle','off');
        surf(Time,Freq, spec)
        view (2)
        shading('interp')
        axis('tight')
        colorbar('eastoutside')
        colormap('parula')
        ylabel('Frequency(Hz)');
        ylim([0 15]);
        caxis([c1 c2])
    catch
        %Error assessment
        if isempty(time) | isempty (spec);set(text1, 'String','Error. Spectrogram matrix is empty.');drawnow;return;end
        if size(Time)~=size(Freq) & size(Time)~=size(spec) & size(Freq)~=size(spec)
        set(text1, 'String','Error. Incoherent dimensions of the Spectrogram.');drawnow;return;
        end    
    end
end

%==========================================================================
%%%%Save the rms plot at different frequency analysis%%%%%%%%%%%%%%%%%%%%%%
function save2(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    try
        ii=slidval;
        amp=cell2mat(AMPP);
        time=datetime(cell2mat(TIME_AMP),'ConvertFrom','datenum');%Convert in datetime format
        amp=amp(:,ii);%Select the ii-th block of the rms
        time=time(:,ii);%Select the ii-th block of the time
        rms_mean=movmean(amp,hours(factor),'SamplePoints',time);%Moving average of the rms according the factor value
        %Plot the normalized amplitude of the rms
        new_figure_2=figure('Name','RMS plot','NumberTitle','off');
        plot(time,rms_mean./max(rms_mean))
        ylabel('Rms (a.u.)');
        grid on 
        grid minor
        legend([num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
        ylim([0 1])
    catch
        %Error assessment
        if isempty(amp) | isempty (time);set(text1, 'String','Error. Rms matrix is empty.');drawnow;return;end
        if size(amp,1)~=size(time,1);set(text1, 'String','Error. Incoherent dimensions of the Rms array.');drawnow;return;end
    end
end

%==========================================================================
%%%%Save the results of the spectrogram%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savespec(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Validation of the results
    if isempty(params.save);set(text1, 'String', sprintf('Error. No output directory.'));drawnow;return;end  
    if isempty(cell2mat(SSTOT));set(text1, 'String', sprintf('Error. Spectrogram matrix is empty.'));drawnow;return;end 
    if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
    np=params.np;
    staz=params.sta_ref;
    %Loop through the number of Input files (ff)
    for kk=1:size(SSTOT,2)
        time=TIME_SPEC{kk};spec=SSTOT{kk};
        [Time Freq]=meshgrid(time,(nyq/(params.np/2)):(nyq/(params.np/2)):nyq);
        %Compute the juliand day
        jday=day(datetime(time(1),'ConvertFrom','datenum'),'dayofyear');
        if jday<10;jday=['00' num2str(jday)];elseif jday>=10 & jday<100;jday=['0' num2str(jday)];elseif jday>=100;jday=[num2str(jday)];end
        %Create the Output directory and filename
        pathout=[params.save '/SPEC/' jday '/'];if ~exist(pathout, 'dir');mkdir(pathout);end 
        name=ff{kk};[~,name,~] = fileparts(name);
        output=[pathout name params.outputformat];
        %Create a table of the results
        data=table(Time,Freq,spec);data.Properties.VariableNames={'Time','Freq', 'Spectrum'};
        data.Properties.VariableUnits = {'datenum','Hz',' '};
        %Save the results
         switch params.outputformat
                case'.mat'
                    save (output,'data');
                case'.txt'
                    writetable(data,output)
          end
    end
    %Command message
    set(text1, 'String', sprintf(['SPEC' params.outputformat '-Observed data are saved.'])); 
    drawnow
end

%==========================================================================
%%%%Save the results of the rms%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saverms(~,~) 
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Validation of the results
    if isempty(params.save);set(text1, 'String', sprintf('Error. No output directory.'));drawnow;return;end  
    if isempty(cell2mat(AMPP));set(text1, 'String', sprintf('Error. Rms matrix is empty.'));drawnow;return;end
    if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
    step=step_rms;
    staz=params.sta_ref;
    %Loop through the frequency analysis
    for jj=1:size(AMPP,2)
        %Loop through the number of Input files (ff)
        for kk=1:size(AMPP,1)
            time=TIME_AMP{kk,jj};amp=AMPP{kk,jj};
            %Compute the juliand day
            jday=day(datetime(time(1),'ConvertFrom','datenum'),'dayofyear');
            if jday<10;jday=['00' num2str(jday)];elseif jday>=10 & jday<100;jday=['0' num2str(jday)];elseif jday>=100;jday=[num2str(jday)];end
            %Create the Output directory and filename
            freq1=F(jj);freq2=F(jj+1);Ffolder=[num2str(freq1) '-' num2str(freq2) 'Hz'];
            name=ff{kk};[~,name,~] = fileparts(name);
            pathout=[params.save '/RMS/' Ffolder '/' jday '/'];if ~exist(pathout, 'dir');mkdir(pathout);end
            output=[pathout name params.outputformat];
            %Create a table of the results
            data=table(time,amp);data.Properties.VariableNames={'Time', 'Rms'};
            data.Properties.VariableUnits = {'datenum','a.u.'};
            %Save the results
            switch params.outputformat
                case'.mat'
                    save (output,'data');
                case'.txt'
                    writetable(data,output)
            end
        end
    end
    %Command message
    set(text1, 'String', sprintf(['RMS' params.outputformat '-Observed data are saved.']));
    drawnow
end
end