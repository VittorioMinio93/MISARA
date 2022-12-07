%%%%%%%%%%%%%%%%%%%%% Salped detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params: user parameters
% avfig: Salped figure
% panel_color: colors of the main panel 
% fig_color: colors of the figure
% pan1: control right side panel 
% text1: Window message
% axes1,axes3: axes of the current figure
% slidval: value of dynamic plots 
% dispval: displacement visualization value
% h: UI control parameters
%%%%%%%%%%%h.variables/params.variables/variables%%%%%%%%%%%%%%%%%%%%%%%%%%
% rad,Rad: radio button within a button group 
% rad1, rad2, rad3, rad4,Rad1, Rad2: button group 
%-------------------------------- General settings
%--- General info
% sta_ref: the name of the reference station 
% component, comp: the name of the component/channel
% comp_chan: the station system used (component/channel)
% station: station names
% coord: the path of the coordinate file (.mat) 
% cooord: coordinates of the station of reference
% coord_sys: the coordinate system (metrics/degrees)
% map: the path of the DEM file (.tif)
% X,Y,A: georeferenced DEM grid NxM (UTM system/WGS84 system)
%--- Data file
% data: the Input folder
% save: the Output folder
% outputformat: the format (.mat/.txt) of the output files
% search_sac, date_search: term of research file 
% pathname/pathnam/pathaname2, fname: Input folder and list of filenames
% ff: list of files with their path
% mytrace: seismic trace structure
% dt: sampling interval
% fs: sampling rate
% nyq: nyquist number
% z,z3: amplitude vector of the seismic trace
% nt: number of samples
% startTime, endTime: temporal limits of the traces 
%--- Detection setting
% AMod: type of saving mode (Automatic/Manual)
% wshort: the short-term window or the pre-trigger window (s)
% wlong: the long-term window or the post-trigger window (s)
% f1_det: the minimum frequency for the central band of interest(Hz)
% f2_det: the maximum frequency for the central band of interest(Hz)
% f1_down: the minimum frequency for the lower band of no interest(Hz)
% f2_down: the maximum frequency for the lower band of no interest(Hz)
% f1_up: the minimum frequency for the upper band of no interest(Hz)
% f2_up: the maximum frequency for the upper band of no interest(Hz)
% treshold: the threshold value of the detection
% lwin: spectrogram window (s)
% nover:  spectrogram window overlap (s)
% fpoint: frequency samples of spectrogram 
% wpol: polarization/rms window (s)
% nratio: number of triggers detected
% time_trigger,TIME_DET: trigger times (s)
% TP_sampl: trigger times (samples)
% juliand, juliandnew: julian day
% CF,VAL: Characteristic function
% TP, TPP, Triggers, TIME_sel: selelction of trigger times (s)
% NEWTIME_SEL: corrected trigger times (s) 
% Name,NAME, Name_sel: file name vector/matrix
% family, Family: Family label 
% Z,Y,YY: extracted waveform vectors
% Z,N,E: the three components
% amp: amplitude of the rms 
% ff1: frequencies of spectrums
% S: spectrums
% az: azimuth vectors
% incid: incidence vectors
% rect: rectilinearity vectors
% sis: three components matrix 
% X_tot, Y_tot, Z_tot: Particle motions matrixs

function salpedDetection()
%% Declaration global parameters
clc
global params
comp=params.component;
component=[];
dispval=0;AMod='Aut';
ff= cell(1);
mytrace=[];
VAL=[];
Family='F1';
Triggers=[];
TIME_DET=[];
NAME=[];NEWTIME_SEL=[];
slidval=1;
dt=[];
fs=[];
YY=[];
Namesel=[];
pathname=[];pathnam=[];pathname2=[];wlong=[];wshort=[];f1=[];f2=[];
%% Salped figure
fig_color= [.94 .94 .94];
panel_color= [.97 .97 .97];

avfig= figure('Name', 'Salped analysis', 'position',get(0,'screensize'), ...
    'NumberTitle','off','units','normalized', 'color', fig_color,'MenuBar','none');

%Window message to user
text1= uicontrol('Style','edit','FontName','Arial','Fontsize',10,'Fontweight','bold',...
                'String','Window Message','Max',5,'Min',0,'HorizontalAlignment','left',...
                'BackgroundColor', 'w','Enable','inactive',...
                'units','normalized','Position',[.05 .02 .72 .05]);

%% Axeses
% Triggers axes
axes1= axes('position',[.05 .56 .72 .40]);
hold(axes1,'on')
     grid(axes1,'on')
     grid(axes1,'minor')
      ylabel(axes1,'CF')
      
% Waveforms axes
axes3= axes('position',[.05 .14 .72 .36]);
hold(axes3,'on')
     grid(axes3,'on')
     grid(axes3,'minor')
     ylabel(axes3,'Amplitude')

%% Control right side panel       
pan1= uipanel(avfig,'visible','on','Position',[.81 .03 .185 .96],...
    'BackgroundColor', panel_color);

uipanel(avfig,'visible','on','Position',[.81 .45 .185 .01],...
    'BackgroundColor', fig_color);

%% Text boxes
%-------------------------------- Salped setting 
uicontrol(pan1,'Style','text', 'String','Salped setting',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.55 .77 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Win1 (s)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .74 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Win2 (s)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .69 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','F1 (Hz)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .64 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','F2 (Hz)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .57 .2 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Thr',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .52 .2 .04],...
    'BackgroundColor',panel_color);

%-------------------------------- SPEC/Pol setting 
uicontrol(pan1,'Style','text', 'String','SPEC/Pol setting',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.03 .77 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Wspe (s)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .74 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Over (s)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .69 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Fpoint',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .62 .2 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Wpol (s)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .57 .2 .04],...
    'BackgroundColor',panel_color);
%-------------------------------- Plot setting
uicontrol(pan1,'Style','text', 'String','Plot setting',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.03 .50 .4 .04],...
    'BackgroundColor',panel_color);

%% Editable text 
%-------------------------------- Salped setting
h.wshort= uicontrol(pan1,'Style','edit', 'String',params.wshort,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .73 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.wlong= uicontrol(pan1,'Style','edit', 'String',params.wlong,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .68 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.f1= uicontrol(pan1,'Style','edit', 'String',params.f1_det,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .63 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.f2= uicontrol(pan1,'Style','edit', 'String',params.f2_det,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .58 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.thr= uicontrol(pan1,'Style','edit', 'String',params.treshold,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .53 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);

%-------------------------------- SPEC/Pol setting 
h.lwin= uicontrol(pan1,'Style','edit', 'String','1.28',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.23 .73 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.nover= uicontrol(pan1,'Style','edit', 'String','1.20',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.23 .68 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.fpoint= uicontrol(pan1,'Style','edit', 'String','512',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.23 .63 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.wpol= uicontrol(pan1,'Style','edit', 'String','2',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.23 .58 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);

h.family= uicontrol(pan1,'Style','popup', ...
    'String',{'F1','F2','F3','F4','No Fam'},'Value',1,...
    'HorizontalAlignment','left','Enable','on',...
    'Units','normalized','Position',[.68 .40 .3 .1],...
    'BackgroundColor','w','callback',@setfamily);
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
    'bordertype','none','Position',[.08 .82 .7 .05]);
set(h.rad,'SelectionChangeFcn',@radcbk);

h.rad1 = uicontrol( h.rad, 'Style','Radio','String','Z',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.0 0 .16 1],'HandleVisibility','off');
h.rad2 = uicontrol( h.rad, 'Style','Radio','String','N',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.18 0 .16 1],'HandleVisibility','off');
h.rad3 = uicontrol( h.rad, 'Style','Radio','String','E',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.36 0 .16 1],'HandleVisibility','off');
h.rad4 = uicontrol( h.rad, 'Style','Radio','String','F',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.54 0 .16 1],'HandleVisibility','off');

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

% Editable text
uicontrol(pan1,'Style','text', 'String','Save Mode',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.63 .85 .4 .05],...
    'BackgroundColor',panel_color);
% Save Mode selection radio buttons
h.Rad= uibuttongroup(pan1,'units','normalized','BackgroundColor',panel_color,...
    'bordertype','none','Position',[.63 .82 .9 .05]);
set(h.Rad,'SelectionChangeFcn',@rad2cbk);
h.Rad1 = uicontrol( h.Rad, 'Style','Radio','String','Aut',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.0 0 .2 1],'HandleVisibility','off');
h.Rad2 = uicontrol( h.Rad, 'Style','Radio','String','Man',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.21 0 .2 1],'HandleVisibility','off');

%Set the save mode (Automatic /Manual)
if strcmp(AMod,'Aut')
    set (h.Rad1,'value',1);
    set (h.Rad2,'value',0);

elseif strcmp(AMod,'Man ')
    set (h.Rad1,'value',0);
    set (h.Rad2,'value',1);
end

%% Buttons
% Open files button
uicontrol('style','pushbutton', 'string','Open files','units','normalized',...
    'position',[.815 .9 .08 .0494], 'callback',@openfnc);

% Open folder button
uicontrol('style','pushbutton', 'string','Open folder','units','normalized',...
    'position',[.91 .9 .08 .0494], 'callback',@openfolder);

% Logarithmic scale checkbox
uicontrol('Style','checkbox','BackgroundColor', panel_color,...
                'String','Log Axis','units','normalized','Tag','logcb',...
                'Value',0,'Position',[.82 .48 .05 .03],'callback',@logscale);
% Displacement visualization checkbox
uicontrol('Style','checkbox','BackgroundColor', panel_color,...
                'String','Disp/Vel','units','normalized','Tag','viscb',...
                'Value',0,'Position',[.88 .48 .05 .03],'callback',@dispvisuale);
            
% Editable texts
uicontrol(pan1,'Style','text', 'String','View all triggers',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.03 .35 .4 .07],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Select triggers',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.55 .35 .4 .07],...
    'BackgroundColor',panel_color);

% Calculation and save button
uicontrol('style','pushbutton', 'string','Calculation','units','normalized',...
    'position',[.815 .32 .08 .0494], 'callback',@calc); 
uicontrol('style','pushbutton', 'string','Calculation','units','normalized',...
    'position',[.91 .32 .08 .0494], 'callback',@calc2);
uicontrol('style','pushbutton', 'string','Save Plot','units','normalized',...
    'position',[.815 .25 .08 .0494], 'callback',@save1); 
uicontrol('style','pushbutton', 'string','Save Plot','units','normalized',...
    'position',[.91 .11 .08 .0494], 'callback',@save2);
uicontrol('style','pushbutton', 'string','Plot Spec','units','normalized',...
    'position',[.91 .25 .08 .0494], 'callback',@save3);
uicontrol('style','pushbutton', 'string','Save Data','units','normalized',...
    'position',[.815 .18 .08 .0494], 'callback',@savetrig);
uicontrol('style','pushbutton', 'string','Plot Pol','units','normalized',...
    'position',[.91 .18 .08 .0494], 'callback',@save4);
uicontrol('style','pushbutton', 'string','Save Data','units','normalized',...
    'position',[.91 .04 .08 .0494], 'callback',@savesel);

% Next and pre button
uicontrol('style','pushbutton', 'string','Next','units','normalized',...
    'position',[.776 .4 .03 .04], 'callback',@stlt2);
uicontrol('style','pushbutton', 'string','Pre','units','normalized',...
    'position',[.776 .46 .03 .04], 'callback',@stlt1);

%% Other functions (nested)
%%%%Return the list of files manually selected%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function openfnc(~,~)
    VAL=[];Triggers=[];TIME_DET=[];NAME=[];NEWTIME_SEL=[];YY=[];Namesel=[];%Initialize the main vectors
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);cla(axes3);%Clear the axes
    tt= 0;pathnam=[];pathname2=[];
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
    VAL=[];Triggers=[];TIME_DET=[];NAME=[];NEWTIME_SEL=[];YY=[];Namesel=[];%Initialize the main vectors
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);cla(axes3);%Clear the axes
    pathname=[];
    pathnam= uigetdir(params.data);
    pathnam=[pathnam '/'];
    search_sac=strcat(params.sta_ref,'.BH',comp,'*.mat');
    component=comp;
    if pathnam~=0
    % Get all files in the directory and subdirectories
    [ff,pathname2]= getAllFiles(pathnam,search_sac);
    ff=ff'; pathname2=pathname2';
    
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
        drawnow;
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
%%%%%%Return the type of saving method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rad2cbk(~,eventdata)
    
    AMod = get(eventdata.NewValue,'String');   
    
end

%==========================================================================
%%%%Calculation and plot of the trigger times%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calc(~,~)
    YY=[];%Initialize the extracted waveforms vector
    NEWTIME_SEL=[];%Initialize the corrected trigger time vector
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);cla(axes3);%Clear the axes
    try
        %Acquisition and validation of the analysis parameters
        f1=str2num(get(h.f1,'string'));
        f2=str2num(get(h.f2,'string'));
        if f1<=0 | f2<=f1;;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        if f2<=0 | f2<=f1;;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        treshold=str2num(get(h.thr,'string'));
        %Command message
        textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
        set(text1, 'String', textLabel); 
        drawnow
        %Reading .mat files
        [~, len]=size (ff);%Number of files selected
        TIME_DET=[];%Initialize the detection time vector
        NAME=[];%Initialize the file names vector
        VAL=[];%Initialize the characteristic function vector
        Triggers=[];%Initialize the trigger times vector
        %Loop through the number of files selected
        for i=1:len
            load(ff{i});%load the i-th file  
            z=mytrace.data;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate; nyq=0.5/dt; nt=mytrace.sampleCount;
            startTime=mytrace.startTime;endTime=mytrace.endTime;
            %Rescale and remove the best straight-fit line
            z=z-mean(z);
            z=detrend(z);
            %Apply the Salped algorithm
            [CF,tCF]=SALPEDfunction(z,fs,f1,f2,params.f1_down,params.f2_down,params.f1_up,params.f2_up);
            % Selection of the trigger times and filenames over the threshold
            [CFmax,locs]=findpeaks(CF,'MinPeakDistance',round(fs),'MinPeakHeight',treshold);
            %Calculation of the trigger times in datetime format
            timemain=datetime(startTime,'ConvertFrom','datenum');%Convert in datetime format
            tCF2=tCF;
            Time=timemain+seconds(tCF2);
            if isempty(CFmax)
                TP=NaN;
            else
                TP=tCF2(locs);
                name=ff{i};%Select the i-th file name 
                %Loop through the number of triggers
                for bb=1:length(TP)
                    TPP=timemain+seconds(TP(bb));TPP.Format='dd-MMM-yyyy HH:mm:ss.SSS';
                    %Concatenate the results of the trigger selection
                    NAME=[NAME;name];
                    Triggers=[Triggers;TPP];
                end
            end
            %Concatenate the results of Salped
            VAL=[VAL;CF];
            Time.Format='dd-MMM-yyyy HH:mm:ss.SSS';
            TIME_DET=[TIME_DET;Time'];
        end
        %Plottinng of the trigger times
        ploting1(TIME_DET,VAL,treshold,axes1)
        %Set the x-axis limits
        xlim(axes1,[min(TIME_DET) max(TIME_DET)])
        %Command message
        set(text1, 'String','Calculation finished.');
        drawnow
    catch
        %Error assessment 
        if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
        if isempty(z);set(text1, 'String','Error. No correct signals data.');drawnow;return;end
        if f1>=fs/2 | f2>fs/2 | f2<=f1;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        if params.f1_down<=0 | params.f2_down<=0 | params.f1_up<=0 | params.f2_up<=0 ;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        if params.f1_down>=fs/2 | params.f2_down>fs/2 | params.f2_down<=params.f1_down;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        if params.f1_up>=fs/2 | params.f2_up>fs/2 | params.f2_up<=params.f1_up;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    end 
end

%==========================================================================
%%%%%%%%%%%%%%%Extraction of the waveforms%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calc2(~,~)
    try
        text1.String=[];set(text1);drawnow;%Clear the Window message
        cla(axes3)%Clear the axes
        %Acquisition and validation of the analysis parameters
        f1=str2num(get(h.f1,'string'));
        f2=str2num(get(h.f2,'string'));
        if f1<=0 | f2<=f1;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        if f2<=0 | f2<=f1;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        if f1>=fs/2 | f2>fs/2;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        wshort=str2num(get(h.wshort,'string'));if wshort<=0;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
        wlong=str2num(get(h.wlong,'string'));if wlong<=0;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
        %Command message
        textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
        set(text1, 'String', textLabel); 
        drawnow
        %Reading the selected trigger times and filenames
        Timesel=Triggers;
        Namesel=NAME;
        %Control and set the Input folder
        if ~isempty(pathnam);pathname=pathname2{1};pathname=[pathname '/'];end;path=pathname(1:end-4);
        ntrig=length(Triggers);%number of triggers selected
        YY=cell(1,ntrig);%Initialize the extracted waveforms vector
        NEWTIME_SEL=cell(ntrig,1);%Initialize the corrected trigger time vector
        %Loop through the number of selected triggers
        for kk=1:ntrig
            %Calculation of the trigger time in samples
            TP_sampl=(minute(Timesel(kk))*60+round(second(Timesel(kk))))*round(fs);
            %Select and load the kk-th files
            file_sac=Namesel(kk,:);
            load(file_sac); 
            z=mytrace.data;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate; nyq=0.5/dt; nt=mytrace.sampleCount;
            startTime=mytrace.startTime;endTime=mytrace.endTime;
            %Rescale and remove the best straight-fit line
            z=z-mean(z);
            z=detrend(z);
            w1=round(wshort*fs);
            w2=round(wlong*fs);Z=[];
            %Extracting the waveforms on the basis of the trigger time  and signal files available
            if TP_sampl-w1>0 & TP_sampl+w2<=length(z)
                %Extraction of the waveform
                Z=z(TP_sampl-w1:TP_sampl+w2);
            elseif TP_sampl-w1<0 & TP_sampl+w2<=length(z)
                %Search the new time of reference
                newtime=Timesel(kk)-1/24;
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
                newfile=strcat(path,'/',juliandnew,'/',params.sta_ref,'.BH',component,'.',num2str(year(newtime)),'.',juliandnew,'.',newtime2,'0000.mat');
                %Control the directory of the new signal file 
                p=dir(newfile);
                if ~isempty(p)
                    %Load the trace
                    load(newfile); 
                    z3=mytrace.data;nt=mytrace.sampleCount;
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
                newtime=Timesel(kk)+1/24;
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
                newfile=strcat(path,'/',juliandnew,'/',params.sta_ref,'.BH',component,'.',num2str(year(newtime)),'.',juliandnew,'.',newtime2,'0000.mat');
                %Control the directory of the new signal file 
                p=dir(newfile);
                if ~isempty(p)
                    %Load the trace
                    load(newfile); 
                    z3=mytrace.data;nt=mytrace.sampleCount;
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
            %Filter the seismic traces with a band-pass
            [Bf,Af]=butter(2,[f1 f2]/(fs/2));
            Y=filter(Bf,Af,Z);
            %Rescale and remove the best straight-fit line
            Y=Y-mean(Y);
            Y=detrend(Y);
            %Selection of the new trigger times (time index of the maximum amplitude)
            [~,imax]=max(Y);
            T_sampl=Timesel(kk);
            t1=T_sampl-seconds(wshort);
            t2=T_sampl+seconds(wlong);
            temp=[t1:seconds(dt):t2];
            [~,idx]=max(Y);
            %Concatenate the results
            newtrig=temp(idx);
            NEWTIME_SEL(kk)= {datenum(newtrig)};
            YY(kk)={Y};
        end
        %Set the block of results
        slidval=ntrig;
        %Error control
        if isempty(cell2mat(YY));set(text1, 'String','Error. Signals matrix is empty.');drawnow;return;end
        %Plotting of the extracted waveforms
        ploting2(slidval,axes3)
        %Command message
        set(text1, 'String','Calculation finished.');
        drawnow
    catch
        %Error assessment
        if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
        if isempty(Triggers);set(text1, 'String','Error.Triggers matrix is empty.');drawnow;return;end
        if isempty(cell2mat(YY));set(text1, 'String','Error. Signals array is empty.');drawnow;return;end
        if isempty(Z);set(text1, 'String','Error. No correct signals data.');drawnow;return;end     
    end    
end

%==========================================================================
%%%%%%%%%%%%%%%%%Plotting of the trigger times%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting1(time,val,treshold,ax)
    try
        cla(ax)%Clear the axes
        hold(ax,'on')
        %Plot the characteristic function values
        plot(ax,time,val)
        %Plot the threshold line
        plot(ax,time,ones(length(val),1).*treshold)
        xlim(ax,[min(time) max(time)])
        ylabel(ax,'CF')
        hold(ax,'off')
    catch
        %Error assessment 
        if size(val,1)~=size(time,1);set(text1, 'String','Error. Incoherent dimensions of the Triggers array.');drawnow;return;end
    end   
end

%==========================================================================
%%%%%%%%%%%%Plotting of the extracted waveforms%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting2(ii,ax)
    try
        time=[];Z=[];%Initialize the waveform and the time vector
        cla(ax)%Clear the axes
        %%Select the ii-th trigger time and set temporal limits of the extracted waveform 
        TP_sampl=Triggers(ii);
        t1=TP_sampl-seconds(wshort);
        t2=TP_sampl+seconds(wlong);
        Y=cell2mat(YY(:,ii));%Select the ii-th waveform 
        %Set the waveform in displacement/velocity
        if dispval==1
            %Expressing  in displacement and normalizing 
            Z=cumtrapz(Y);
            Z=Z./max(Z);
        elseif dispval==0
            %Expressing  in velocity and normalizing
            Z=Y;
            Z=Z./max(Z);       
        end
        %Plot of the extracted waveform 
        time=[t1:seconds(dt):t2];
        plot(ax,time,Z)
        ylabel(ax,'Amplitude');
        t_sel=datetime(cell2mat(NEWTIME_SEL(ii)),'ConvertFrom','datenum');%Convert in datetime format
        legend(ax,datestr(t_sel,'dd-mmm-yyyy HH:MM:SS.FFF'))
        xlim(ax,[t1,t2])
        ylim(ax,[min(Z) max(Z)])
    catch
        %Error assessment 
        if size(Z,1)~=size(time,1);set(text1, 'String','Error. Incoherent dimensions of the Signals array.');drawnow;return;end
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
%%%%%%%%%%%%%%%%%%%%%Set the visualization mode %%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispvisuale(h, ~)
    if ~isempty(cell2mat(YY))
        if (get(h,'Value') == get(h,'Max'))
            %Command message
            set(text1, 'String', 'Plotting in displacement.');
            %Set the visualization mode value
            dispval=1;
            % Plot the seismic traces in displacement
            ploting2(slidval,axes3)
        else
            %Command message
            set(text1, 'String', 'Plotting in velocity.');
            %Set the visualization mode value
            dispval=0;
            % Plot the seismic traces in velocity
            ploting2(slidval,axes3)
        end
    else
        if (get(h,'Value') == get(h,'Max'))
            %Error message
            set(text1, 'String', 'Error. Signals matrix is empty.');
            drawnow
            %Set the visualization mode value
            dispval=1;
            return
        else
            %Error message
            set(text1, 'String', 'Error. Signals matrix is empty.');
            drawnow
            %Set the visualization mode value
            dispval=0;
            return
        end
    end
end

%==========================================================================
%%%%Refresh the waveform axes through the pre button%%%%%%%%%%%%%%%%%%%%%%%
function stlt1(~,~)
   %Error control
   if isempty(Triggers);set(text1, 'String', 'Error. Signals matrix is empty.');return;end
   %Decrease the slidval value
   slidval=slidval-1;if slidval<=0;slidval=1;end
   %Plot the selected waveform
   ploting2(slidval,axes3);
        
end

%==========================================================================
%%%%Refresh the waveform axes through the next button%%%%%%%%%%%%%%%%%%%%%%
function stlt2(~,~)
   %Error control
   if isempty(Triggers);set(text1, 'String', 'Error. Signals matrix is empty.');return;end
   %Increase the slidval value
   slidval=slidval+1;if slidval>length(Triggers);slidval=length(Triggers);end
   %Plot the selected waveform 
   ploting2(slidval,axes3);        
end

%==========================================================================
%%%%Control and reset automatically analysis parameters%%%%%%%%%%%%%%%%%%%%
function deflimits(~,~)
    %Acquisition analysis parameters
    f1=str2num(get(h.f1,'string'));
    f2=str2num(get(h.f2,'string'));
    wshort=str2num(get(h.wshort,'string'));
    wlong=str2num(get(h.wlong,'string'));
    treshold=str2num(get(h.thr,'string'));
    lwin=str2num(get(h.lwin,'string'));
    nover=str2num(get(h.nover,'string'));
    fpoint=str2num(get(h.fpoint,'string'));
    wpol=str2num(get(h.wpol,'string'));
    
    %Control and eventually reset parameters
    if (size(wshort,1)==0);set(h.wshort,'string',params.wshort);end
    if (size(wlong,1)==0);set(h.wlong,'string',params.wlong);end
    if (size(f1,1)==0);set(h.f1,'string',params.f1_det);end
    if (size(f2,1)==0);set(h.f2,'string',params.f2_det);end
    if (size(treshold,1)==0);set(h.thr,'string',params.treshold);end
    if (size(lwin,1)==0);set(h.lwin,'string','1.28');end
    if (size(nover,1)==0);set(h.nover,'string','1.20');end
    if (size(fpoint,1)==0);set(h.fpoint,'string','512');end
    if (size(wpol,1)==0);set(h.wpol,'string','2');end
end

%==========================================================================
%%%%%%%%%%%%%%%%%%%%%Set the family label %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setfamily(h, ~)
    if get(h,'Value')==1
        Family='F1';
    elseif get(h,'Value')==2
        Family='F2';
    elseif get(h,'Value')==3
        Family='F3';
    elseif get(h,'Value')==4
        Family='F4';
    elseif get(h,'Value')==5
        Family='No Fam';    
    end
    
end

%==========================================================================
%%%%%%%%%%%%%%%%%%%Save triggers plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save1(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Error control
    if isempty(TIME_DET);set(text1, 'String', 'Error. Triggers matrix is empty.');return;end
    try
        %Plot the characteristic function values
        new_figure=figure('Name','View all triggers','NumberTitle','off');
        hold on
        plot(TIME_DET,VAL)
        %Plot the threshold line
        plot(TIME_DET,ones(length(VAL),1).*treshold)
        xlim([min(TIME_DET) max(TIME_DET)])
        hold off
        ylabel('CF')
        grid on;grid minor;
    catch
        %Error assessment
        if size(TIME_DET,1)~=size(VAL,1);set(text1, 'String','Error. Incoherent dimensions of the Triggers array.');drawnow;return;end
    end  
end

%==========================================================================
%%%%%%%%%%%%%%%%%%%Save the waveform  plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save2(~,~)
    Z=[];time=[];%Initialize the waveform and the time vector
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Acquisition and validation of the analysis parameters
    wshort=str2num(get(h.wshort,'string'));if wshort<=0;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
    wlong=str2num(get(h.wlong,'string'));if wlong<=0;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
    if f1<=0 | f2<=f1;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    if f2<=0 | f2<=f1;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    if f1>=fs/2 | f2>fs/2;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    try 
        %%Select the ii-th trigger time and set temporal limits of the extracted waveform 
        ii=slidval;
        t1=Triggers(ii)-seconds(wshort);
        t2=Triggers(ii)+seconds(wlong);
        Y=cell2mat(YY(:,ii));%Select the ii-th waveform 
        %Set the waveform in displacement/velocity
        if dispval==1
            %Expressing  in displacement and normalizing 
            Z=cumtrapz(Y);
            Z=Z./max(Z);
        elseif dispval==0
            %Expressing  in velocity and normalizing
            Z=Y;
            Z=Z./max(Z);       
        end
        %Error control
        if isempty(Z);set(text1, 'String','Error. Signals matrix is empty.');drawnow;return;end
        time=[t1:seconds(dt):t2];
        %Plot of the extracted waveform 
        new_figure_2=figure('Name','Select Trigger','NumberTitle','off');
        plot(time,Z)
        ylabel('Amplitude');
        grid on;
        grid 'minor';
        t_sel=cell2mat(NEWTIME_SEL(ii));
        legend(datestr(t_sel,'dd-mmm-yyyy HH:MM:SS.FFF'))
        xlim([t1,t2])
        ylim([min(Z) max(Z)])
    catch
        %Error assessment
        if isempty(Z);set(text1, 'String','Error. Signals matrix is empty.');drawnow;return;end
        if size(Z,1)~=size(time,1);set(text1, 'String','Error. Incoherent dimensions of the Triggers array.');drawnow;return;end    
    end
end

%==========================================================================
%%%%%%%%%%%%%%%Calculation and plotting of the spectrogram%%%%%%%%%%%%%%%%%  
function save3(~,~)
    Z=[];time=[];%Initialize the waveform and the time vector
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Acquisition and validation of the analysis parameters
    wshort=str2num(get(h.wshort,'string'));if wshort<=0;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
    wlong=str2num(get(h.wlong,'string'));if wlong<=0;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
    if f1<=0 | f2<=f1;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    if f2<=0 | f2<=f1;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    if f1>=fs/2 | f2>fs/2;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    lwin=round(str2num(get(h.lwin,'string'))*fs); if lwin<=0 | lwin>=round((wlong+wshort+dt)*fs) ;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
    nover=round(str2num(get(h.nover,'string'))*fs);if nover<0 | nover>=lwin;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
    fpoint=str2num(get(h.fpoint,'string'));if  fpoint<=0 | fpoint>=100000 ;set(text1, 'String','Error. Invalid analysis parameters.');drawnow;return;end
    try 
        %%Select the ii-th trigger time and set temporal limits of the extracted waveform 
        ii=slidval;
        t1=Triggers(ii)-seconds(wshort);
        t2=Triggers(ii)+seconds(wlong);
        Y=cell2mat(YY(:,ii));%Select the ii-th waveform 
        %Calculation of the spectrogram 
        [S, ff1, ts] = spectrogram(Y,lwin,nover,fpoint,fs,'yaxis');
        %Set the temporal limits of the spectrogram
        tt1=t1+seconds(ts(1));tt2=t1+seconds(ts(end));
        %Normalization of the spectrum
        S=abs(S);S = S/max(S(:))+1e-5; 
        time=[t1:seconds(dt):t2];
        new_figure_2=figure('Name','Spectrogram vs. Signal','NumberTitle','off');
        %Plot the spectrogram 
        subplot(2,1,1)
        surf(ts,ff1,S);
        shading interp
        colormap parula
        cb=colorbar('eastoutside');
        cb.Position = cb.Position + [.1 0 0 0];
        xlim([min(ts),max(ts)])
        set(gca,'XTickLabel',[]);
        ylabel('Frequency (Hz)')
        view([0 90])
        ylim([0 5])
        subplot(2,1,2)
        %Set the waveform in displacement/velocity
        if dispval==1
            %Expressing  in displacement and normalizing 
            Z=cumtrapz(Y);
            Z=Z./max(Z);
        elseif dispval==0
            %Expressing  in velocity and normalizing 
            Z=Y;
            Z=Z./max(Z);       
        end
        %Plot the extracted waveform 
        plot(time,Z)
        ylabel('Amplitude');
        grid on;
        grid 'minor';
        t_sel=cell2mat(NEWTIME_SEL(ii));
        legend(datestr(t_sel,'dd-mmm-yyyy HH:MM:SS.FFF'))
        xlim([tt1,tt2]) 
        ylim([min(Z) max(Z)])
    catch
        %Error assessment 
        if isempty(Z);set(text1, 'String','Error. Signals matrix is empty.');drawnow;return;end
        if size(Z,1)~=size(time,1);set(text1, 'String','Error. Incoherent dimensions of the Triggers array.');drawnow;return;end   
    end
end

%==========================================================================
%%%%%%%%%%%%%%%Calculation and plotting of the particle motions%%%%%%%%%%%%
function save4(~,~) 
    %Acquisition and validation of the analysis parameters
    wshort=str2num(get(h.wshort,'string'));if wshort<=0;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
    wlong=str2num(get(h.wlong,'string'));if wlong<=0;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
    if f1<=0 | f2<=f1;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    if f2<=0 | f2<=f1;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    if f1>=fs/2 | f2>fs/2;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    wls=str2num(get(h.wpol,'string'));if wls<=0 | wls>=(wlong+wshort+dt);set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
    try
        A=[];X=[];Y=[];xLimits=[];yLimits=[];%Initialize the DEM variables  
        cooord=[]; utm_zone=[];%Initialize the station coordinates variables   
        text1.String=[];set(text1);drawnow;%Clear the Window message
        %Execute one of two groups of routines (Channel/Component system)
        switch params.comp_chan
            case'Comp'
                %Select the ii-th extracted waveform 
                ii=slidval;
                ZZ=cell2mat(YY(:,ii));               
                %Select the ii-th trigger time 
                time=Triggers(ii);
                %Set the path of the station  coordinates file 
                path=params.coord;
                %Control and set the Input folder
                if ~isempty(pathnam);pathname=pathname2{1};pathname=[pathname '/'];end;mainpath=pathname(1:end-4);
                %Set the temporal limits of waveform 
                t1=time-seconds(wshort);
                t2=time+seconds(wlong);
                %Select the three components of  the waveforms detected
                [Z,N,E]=sselect3comp(time,Namesel(ii,:),mainpath,fs,wlong,wshort,params.sta_ref,text1);
                %Create the three component matrix
                sis(:,1)=Z;sis(:,2)=N;sis(:,3)=E;
                %Application of the polarization analysis
                [az,incid,rect,LL,TL] = filter_taper(sis,fs,f1,f2,wls,wls);
                Time=t1:seconds(wls):t2;
                Time=Time(1:length(az))';
                %Application of the Particle motion algorithm
                [X_tot,Y_tot,Z_tot,cooord]=ParticleMotion(path,sis,fs,f1,f2,wls,params.sta_ref,params.coord_sys,text1);
                %Error control
                if isempty(params.coord) | isempty(cooord);set(text1, 'String','Error. No coordinate file directory or No correct coordinate file.');drawnow;return;end
                %Application of the RMS algorithm
                amp=RMS(Z,fs,dt,f1,f2,wls);
                %Control the DEM file availability
                if ~isempty(params.map)
                    %Execute one of two groups of routines (UTM/wgs84 system)
                    switch params.coord_sys
                        case 'm'
                            %Load the DEM variables
                            load('mapfile_m.mat')
                            %Set the elevation limits
                            zmin=min(min(A));zmax=max(max(A));
                            %Calcultion of the index of maximum rms value
                            [val,imax]=max(amp);
                            %Set the title with the polarization results
                            tit=cellstr(sprintf('Az=%.2f, Inc=%.2f, Ret=%.2f \n',az(imax),incid(imax),rect(imax)));
                            new_figure_2=figure('Name','Polarization attributes','NumberTitle','off');
                            %Plot the XY section
                            subplot(2,2,1)
                            %Plot the DEM 
                            contour3(X,Y,A,(zmin:0.2:zmax),'k')
                            xlabel('X (km)')
                            ylabel('Y (km)')
                            hold on
                            %Plot the stations
                            scatter3(cooord(1),cooord(2),cooord(3),'dk','filled')
                            %Plot of the particle motions with greater energy
                            plot3(X_tot(:,imax),Y_tot(:,imax),Z_tot(:,imax),'r')
                            hold off
                            view([0 90])
                            title(tit)
                            %Plot the YZ section
                            subplot(2,2,2)
                            %Plot the DEM 
                            contour3(X,Y,A,(zmin:0.2:zmax),'k')
                            ylabel('Y (km)')
                            zlabel('Z (km)')
                            hold on
                            %Plot the stations
                            scatter3(cooord(1),cooord(2),cooord(3),'dk','filled')
                            %Plot of the particle motions with greater energy
                            plot3(X_tot(:,imax),Y_tot(:,imax),Z_tot(:,imax),'r')
                            hold off
                            view([90 0])
                            %Plot the XZ section
                            subplot(2,2,3)                            
                            %Plot the DEM 
                            contour3(X,Y,A,(zmin:0.2:zmax),'k')
                            xlabel('X (km)')
                            zlabel('Z (km)')
                            hold on
                            %Plot the stations
                            scatter3(cooord(1),cooord(2),cooord(3),'dk','filled')
                            %Plot of the particle motions with greater energy
                            plot3(X_tot(:,imax),Y_tot(:,imax),Z_tot(:,imax),'r')
                            hold off
                            view([0 0])
                        case 'degrees'
                            %Load the DEM variables
                            load('mapfile_deg.mat')
                            A=A*1000;%Expressed in meters
                            %Set the elevation limits
                            zmin=min(min(A));zmax=max(max(A));
                            %Calcultion of the index of maximum rms value
                            [val,imax]=max(amp);
                            %Set the title with the polarization results
                            tit=cellstr(sprintf('Az=%.2f, Inc=%.2f, Ret=%.2f \n',az(imax),incid(imax),rect(imax)));
                            new_figure_2=figure('Name','Polarization attributes','NumberTitle','off');
                            %Plot the XY section
                            subplot(2,2,1)
                            %Plot the DEM 
                            contour3(X,Y,A,(zmin:200:zmax),'k')
                            xlabel('Long (deg)')
                            ylabel('Lat (deg)')
                            hold on
                            %Plot the stations
                            scatter3(cooord(2),cooord(1),cooord(3),'dk','filled')
                            %Plot of the particle motions with greater energy
                            plot3(X_tot(:,imax),Y_tot(:,imax),Z_tot(:,imax),'r')
                            hold off
                            view([0 90])
                            title(tit)
                            %Plot the YZ section
                            subplot(2,2,2)
                            %Plot the DEM
                            contour3(X,Y,A,(zmin:200:zmax),'k')
                            ylabel('Lat (deg)')
                            zlabel('Elev (m)')
                            hold on
                            %Plot the stations
                            scatter3(cooord(2),cooord(1),cooord(3),'dk','filled')
                            %Plot of the particle motions with greater energy
                            plot3(X_tot(:,imax),Y_tot(:,imax),Z_tot(:,imax),'r')
                            hold off
                            view([90 0])
                            %Plot the XZ section
                            subplot(2,2,3)
                            %Plot the DEM 
                            contour3(X,Y,A,(zmin:200:zmax),'k')
                            xlabel('Long (deg)')
                            zlabel('Elev(m)')
                            hold on
                            %Plot the stations
                            scatter3(cooord(2),cooord(1),cooord(3),'dk','filled')
                            %Plot of the particle motions with greater energy
                            plot3(X_tot(:,imax),Y_tot(:,imax),Z_tot(:,imax),'r')
                            hold off
                            view([0 0])     
                    end
                else
                    %Execute one of two groups of routines (UTM/wgs84 system)
                    switch params.coord_sys
                        case 'm'
                            %Calcultion of the index of maximum rms value
                            [val,imax]=max(amp);
                             %Set the title with the polarization results
                            tit=cellstr(sprintf('Az=%.2f, Inc=%.2f, Ret=%.2f \n',az(imax),incid(imax),rect(imax)));
                            new_figure_2=figure('Name','Polarization attributes','NumberTitle','off');
                            %Plot the XY section
                            subplot(2,2,1)
                            xlabel('X (km)')
                            ylabel('Y (km)')
                            %Plot the stations
                            scatter3(cooord(1),cooord(2),cooord(3),'dk','filled')
                            %Plot of the particle motions with greater energy
                            plot3(X_tot(:,imax),Y_tot(:,imax),Z_tot(:,imax),'r')
                            grid on 
                            grid minor
                            hold off
                            view([0 90])
                            title(tit)
                            %Plot the YZ section
                            subplot(2,2,2)
                            ylabel('Y (km)')
                            zlabel('Z (km)')
                            hold on
                            %Plot the stations
                            scatter3(cooord(1),cooord(2),cooord(3),'dk','filled')
                            %Plot of the particle motions with greater energy
                            plot3(X_tot(:,imax),Y_tot(:,imax),Z_tot(:,imax),'r')
                            grid on 
                            grid minor
                            hold off
                            view([90 0])
                            %Plot the XZ section
                            subplot(2,2,3)
                            xlabel('X (km)')
                            zlabel('Z (km)')
                            hold on
                            %Plot the stations
                            scatter3(cooord(1),cooord(2),cooord(3),'dk','filled')
                            %Plot of the particle motions with greater energy
                            plot3(X_tot(:,imax),Y_tot(:,imax),Z_tot(:,imax),'r')
                            grid on 
                            grid minor
                            hold off
                            view([0 0])
                        case 'degrees'
                            %Calcultion of the index of maximum rms value
                            [val,imax]=max(amp);
                            %Set the title with the polarization results
                            tit=cellstr(sprintf('Az=%.2f, Inc=%.2f, Ret=%.2f \n',az(imax),incid(imax),rect(imax)));
                            new_figure_2=figure('Name','Polarization attributes','NumberTitle','off');
                            %Plot the XY section
                            subplot(2,2,1)
                            xlabel('Long (deg)')
                            ylabel('Lat (deg)')
                            hold on
                            %Plot the stations
                            scatter3(cooord(2),cooord(1),cooord(3),'dk','filled')
                            %Plot of the particle motions with greater energy
                            plot3(X_tot(:,imax),Y_tot(:,imax),Z_tot(:,imax),'r')
                            grid on 
                            grid minor
                            hold off
                            view([0 90])
                            title(tit)
                            %Plot the YZ section
                            subplot(2,2,2)
                            ylabel('Lat (deg)')
                            zlabel('Elev (m)')
                            hold on
                            %Plot the stations
                            scatter3(cooord(2),cooord(1),cooord(3),'dk','filled')
                            %Plot of the particle motions with greater energy
                            plot3(X_tot(:,imax),Y_tot(:,imax),Z_tot(:,imax),'r')
                            grid on 
                            grid minor
                            hold off
                            view([90 0])
                            %Plot the XZ section
                            subplot(2,2,3)
                            xlabel('Long (deg)')
                            zlabel('Elev(m)')
                            hold on
                            %Plot the stations
                            scatter3(cooord(2),cooord(1),cooord(3),'dk','filled')
                            %Plot of the particle motions with greater energy
                            plot3(X_tot(:,imax),Y_tot(:,imax),Z_tot(:,imax),'r')
                            grid on 
                            grid minor
                            hold off
                            view([0 0])     
                    end
                end
            case 'Chan'
                %Error message
                textLabel = sprintf('Error. No three components files are found.');
                set(text1, 'String', textLabel);
                drawnow
                return
        end
    catch
       %Error assessment
       if isempty(Namesel);set(text1, 'String','Error. No files are loaded.');drawnow;return;end
       if isempty(Triggers);set(text1, 'String','Error.Triggers matrix is empty.');drawnow;return;end
       if isempty(params.coord) | isempty(cooord);set(text1, 'String','Error. No coordinate file directory or No correct coordinate file.');drawnow;return;end
       if isempty(cell2mat(YY));set(text1, 'String','Error. Signals array is empty.');drawnow;return;end
       if f1>=fs/2 | f2>fs/2;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    end
end

%==========================================================================
%%%%Calculation of the rms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function amp=RMS(z,fs,dt,f1,f2,wls)
    %Filter the seismic traces with a band-pass
    fs=(1/dt);fn=fs/2;
    W(1,1)=f1/fn; W(1,2)=f2/fn;
    [a b]=butter(2,W);
    z=filter(a,b,z);
    %Set the analysis parameters
    N=length(z);wl=round(wls/dt);
    i1=0; i2=0;
    k=0;
    %Move the analysis window through the seismic traces
    while i2 < N-wl || i2== N-wl
        k=k+1;
        i1=(1+(k-1)*wl); i2=i1+wl-1;%Limits of the analysis window
        y1=z(i1:i2).*tukeywin(wl);%tapered cosine window
        %Apply the rms algorithm
        amp(k)=sqrt(sum((y1.*conj(y1))/wl));
    end
    amp=amp';
end

%==========================================================================
%%%%%%%%%%%%%%%Save the results of the Salped analysis %%%%%%%%%%%%%%%%%%%%
function savetrig(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Validation of the results
    if isempty(params.save);set(text1, 'String', sprintf('Error. No output directory.'));drawnow;return;end  
    if isempty(Triggers);set(text1, 'String', sprintf('Error. Triggers matrix is empty.'));drawnow;return;end
    %Create the Output directory and filename
    if ~exist([params.save '/Salped/'], 'dir');mkdir([params.save '/Salped/']);end
    output=[params.save '/Salped/' 'Time_triggers' params.outputformat];
    %Create a table of the results
    data=table(datenum(Triggers));data.Properties.VariableNames={'Time'};data.Properties.VariableUnits = {'datenum'};
    %Save the results
    if params.outputformat=='.mat'
        save (output,'data');
        %Command message
        set(text1, 'String', sprintf('Triggers.mat-Observed data are saved.'));
        drawnow
    elseif params.outputformat=='.txt'
        writetable(data,output)
        %Command message
        set(text1, 'String', sprintf('Triggers.txt-Observed data are saved.'));
        drawnow
    end   
end

%==========================================================================
%%%%%%%%%%%%%%%Save the extracted waveforms%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savesel(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Validation of the results
    if isempty(params.save);set(text1, 'String', sprintf('Error. No output directory.'));drawnow;return;end  
    if isempty(params.coord);set(text1, 'String', sprintf('Error. No coordinate file directory or No correct coordinate file.'));drawnow;return;end
    if isempty(cell2mat(NEWTIME_SEL));set(text1, 'String', sprintf('Error. Triggers matrix is empty.'));drawnow;return;end
    if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
    ii=slidval;%Set the slidava value
    path=params.coord;%Set the path of the station coordinates file 
    %Control and set the Input folder
    if ~isempty(pathnam);pathname=pathname2{1};pathname=[pathname '/'];end;mainpath=pathname(1:end-4);
    %Create the Output directory
    pathout=[params.save '/Salped/' Family '/'];
    outputformat=params.outputformat;
    if ~exist(pathout, 'dir');mkdir(pathout);end 
    %Execute one of two groups of routines (Automatic/manual mode)
    switch AMod
        case 'Man'
            %Manual saving 
            time=NEWTIME_SEL{ii};%Select the ii-th trigger time 
            %Execute one of two groups of routines (Channel/Component system)
            switch params.comp_chan
                case 'Chan'
                    %Save the ii-th waveforms detected (channel system)
                     ttime1comp(time,Namesel(ii,:),path,mainpath,fs,wlong,wshort,pathout,outputformat,text1,ii,length(NEWTIME_SEL))            
                case 'Comp'
                    %Save the ii-th waveforms detected (component system)
                    ttime3comp(time,Namesel(ii,:),path,mainpath,fs,wlong,wshort,pathout,outputformat,text1,ii,length(NEWTIME_SEL))
            end
        case 'Aut'
            %Automatic saving
            %Execute one of two groups of routines (Channel/Component system)
             switch params.comp_chan
                case 'Chan'
                    %Loop through the number of triggers 
                    for kk=1:length(NEWTIME_SEL)
                        time=NEWTIME_SEL{kk};%Select the kk-th trigger time 
                        %Save the kk-th waveforms detected (channel system)
                        ttime1comp(time,Namesel(kk,:),path,mainpath,fs,wlong,wshort,pathout,outputformat,text1,kk,length(NEWTIME_SEL))                        
                    end                    
                case 'Comp'
                    %Loop through the number of triggers 
                    for kk=1:length(NEWTIME_SEL)
                        time=NEWTIME_SEL{kk};%Select the kk-th trigger time
                        %Save the kk-th waveforms detected (component system)
                        ttime3comp(time,Namesel(kk,:),path,mainpath,fs,wlong,wshort,pathout,outputformat,text1,kk,length(NEWTIME_SEL))                        
                    end
             end
    end
end
end 