%%%%%%%%%%%%% Calculation of the Coheregram vs. Rms %%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params: user parameters
% avfig: Coheregram figure
% panel_color: colors of the main panel 
% fig_color: colors of the figure
% pan1: control right side panel 
% text1: Window message
% axes1,axes3: axes of the current figure
% f1_lim, f2_lim: limit of y-axis on coheregram axes
% t1_lim, t2_lim: limit of x-axis on coheregram and rms axes
% slidval: value of dynamic plots 
% logval: set the scale of the spectral amplitude
% h: UI control parameters
%%%%%%%%%%%h.variables/params.variables/variables%%%%%%%%%%%%%%%%%%%%%%%%%%
% rad: radio button within a button group 
% rad1, rad2, rad3, rad4: button group 
%-------------------------------- General settings
%--- General info
% sta_ref: the name of the reference station 
% component, comp: the name of the component/channel
% comp_chan: the station system used (component/channel)
% station: station names
% coord: the path of the coordinate file (.mat) 
% coord_sys: the coordinate system (metrics/degrees)
%--- Data file
% data: the Input folder
% save: the Output folder
% outputformat: the format (.mat/.txt) of the output files
% search_sac, date_search: term of research file 
% pathname/pathnam/pathaname2, fname: Input folder and list of filenames
% ff: list of files with their path
% ff2: list of seismic traces directories
% mytrace: seismic trace structure
% dt: sampling interval
% fs: sampling rate
% nyq: nyquist number
% z: amplitude vector of the seismic trace
% nt: number of samples
% startTime, endTime: temporal limits of the traces 
%--- Spectral setting
% np: the dimension of the spectrums
% f_high_pass: the frequency for the high-pass filter (Hz)
% Time_coh, TIME_COH: time vectors of the coheregram 
% CXY, SS, SPE_COR: spectrum vectors
% F,FTOT: frequency vector (Hz)
% FF: frequencies vector of the maximum spectral points (Hz)
% Z1: position of the maximum spectral points
%--- RMS setting
% w_rms: the analysis window for the RMS analysis (s)
% f1: the minimum frequency for the RMS analysis (Hz)
% f2: the maximum frequency for the RMS analysis (Hz)
% step_rms: the frequency step for the RMS analysis (Hz)
% factor: the temporal factor for the moving average of the RMS (hours) 
% Time_rms, TIME_RMS, TIME_AMP: time vectors of the rms
% amp, AMP,  AMPP: amplitude vectors of the rms
% FREQ:  range of frequency analysis (Hz) 

function spectral_coh()
%% Declaration global parameters
clc
global params
comp=params.component;
component=[];
ff= cell(1);
mytrace=[];
TIME_AMP=[];
TIME_COH=[];
AMPP=[];
SPE_COR=[];
FTOT=[];station=[];
FF=[];
Z1=[];
pathname=[];pathnam=[];pathname2=[];
fname=[];
slidval=1;
FREQ=[];
f1=[];f2=[];step_rms=[];factor=[];w_rms=[];
f1_lim=[];f2_lim=[];
t1_lim=[];t2_lim=[];
logval=0;
%% Spectral coherence figure
fig_color= [.94 .94 .94];
panel_color= [.97 .97 .97];

avfig= figure('Name', 'Spectral Coherence', 'position',get(0,'screensize'), ...
    'NumberTitle','off','units','normalized', 'color', fig_color,'MenuBar','none');%%Create the spectral coherence figure

%Window message to user
text1= uicontrol('Style','edit','FontName','Arial','Fontsize',10,'Fontweight','bold',...
                'String','Window Message','Max',5,'Min',0,'HorizontalAlignment','left',...
                'BackgroundColor', 'w','Enable','inactive',...
                'units','normalized','Position',[.05 .02 .72 .05]);


%% Axeses
%Coheregram axes
axes1= axes('position',[.05 .56 .72 .40]);
hold(axes1,'on')
ylim(axes1,[0 15])
ylabel(axes1,'Frequency(Hz)');
colormap(axes1,'parula')
colorbar(axes1,'northoutside');
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
    'bordertype','none','Position',[.08 .82 .9 .05]);

set(h.rad,'SelectionChangeFcn',@radcbk);

h.rad1 = uicontrol( h.rad, 'Style','Radio','String','Z',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.0 0 .15 1],'HandleVisibility','off');
h.rad2 = uicontrol( h.rad, 'Style','Radio','String','N',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.2 0 .15 1],'HandleVisibility','off');
h.rad3 = uicontrol( h.rad, 'Style','Radio','String','E',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.4 0 .15 1],'HandleVisibility','off');
h.rad4 = uicontrol( h.rad, 'Style','Radio','String','F',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.6 0 .15 1],'HandleVisibility','off');

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

% Logarithmic scale checkbox
uicontrol('Style','checkbox','BackgroundColor', panel_color,...
                'String','Log Axis','units','normalized','Tag','logcb',...
                'Value',0,'Position',[.82 .49 .05 .03],'callback',@logscale);
% Logarithmic spectrum checkbox
uicontrol('Style','checkbox','BackgroundColor', panel_color,...
                'String','Log Spec','units','normalized','Tag','logcb',...
                'Value',0,'Position',[.82 .46 .05 .03],'callback',@logspec);
% Editable texts
uicontrol(pan1,'Style','text', 'String','Spectral Coherence',...
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
   TIME_AMP=[];TIME_COH=[];AMPP=[];SPE_COR=[];%Initialize the coheregram and rms vectors
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);cla(axes3);%Clear the axes
    tt= 0;pathnam=[];pathname2=[];
    ff=cell(1);%Initialize the list of files
    search_sac=strcat(params.sta_ref,'.BH',comp,'*.mat');
    [fname, pathname] = uigetfile({search_sac;'*.*'},'File Selector',...
        params.data,'MultiSelect','on');%Selection of the files     
    component=comp;
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
    TIME_AMP=[];TIME_COH=[];AMPP=[];SPE_COR=[];%Initialize the coheregram and rms vectors
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
%%%%Calculation and plot of the coheregram%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calc(~,~) 
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1)%Clear the axes
    try
        station=[];%Initialize the stations names vector
        load(params.coord,'station');%Load the station names vector from coordinates file 
        %Validation of the analysis parameters
        if params.np<=0; set(text1, 'String','Error. Invalid Coherogram parameters.');drawnow;return;end
        if params.f_high_pass<=0;set(text1, 'String','Error. Invalid High pass filter frequency.');drawnow;return;end
        if params.n_sta<=0; set(text1, 'String','Error. Invalid Coherogram parameters.');drawnow;return;end  
        %Command message
        textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
        set(text1, 'String', textLabel); 
        drawnow
        %Reading .mat files
        [~, len]=size (ff);%Number of files selected
        TIME_COH=cell(len,1);%Initialize the time vector
        SPE_COR=cell(1,len);%Initialize the spectral coherence vector
        FTOT=cell(1,len);%Initialize the frequency vector
        %Loop through the number of files selected
        for i=1:len
            date_search=char(ff{i});%return the i-th file
            date_search=date_search(end-18:end-4);%Create the term of research
            %Control and set the Input folder
            if ~isempty(pathnam);pathname=pathname2{i};pathname=[pathname '/'];end
            %Search all seismic traces based on the station names
            ff2=searchAllStaz(pathname,station,component,date_search,params.comp_chan);
            nny=size(ff2);
            nny=nny(1,1);%Number of traces selected
            delta=0;
            %Choose the number of traces
            if  nny<params.n_sta;delta=abs(nny-params.n_sta);end
            ff2=ff2(1:params.n_sta-delta);
            nny=size(ff2);
            nny=nny(1,1);%Number of traces selected
            %Control the calculation through the number of stations
            if nny==params.n_sta
                SS=[];%Initialize the spectrums vector
                Z=[];%Initialize the signals vector
                %Loop through the number of stations used
                for ifile=1:nny
                    filez=strcat(pathname,char(ff2(ifile).name));   
                    load(filez);%load the ifile-th seismic trace 
                    z=mytrace.data;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate; nyq=0.5/dt; nt=mytrace.sampleCount;
                    startTime=mytrace.startTime;endTime=mytrace.endTime;
                    %Rescale and remove the best straight-fit line
                    z=z-mean(z);
                    z=detrend(z);
                    %Filter the seismic traces with a high-pass
                    ff1=params.f_high_pass;
                    fn=fs/2;
                    W=ff1/fn;
                    [x y]=butter(2,W,'high');
                    z=filter(x,y,z);
                    %Concatenate the seismic traces through the columns
                    ntt(ifile)=length(z);
                    Z(1:ntt(ifile),ifile)=z;
                end
                %Control and set the length of the seismic traces
                nt=min(ntt); temp=Z; clear Z; Z(1:nt,:)=temp(1:nt,:); clear temp
                %Loop through all independent pairs of stations
                for is1=1:nny-1
                    for is2=is1+1:nny
                        %Apply the Magnitude-squared coherence algorithm
                        [CXY,F] = mscohere(Z(:,is1),Z(:,is2),hanning(params.np),params.np/2,params.np,fs);
                        %Concatenate the spectrums
                        SS=[SS,CXY];
                    end
                end 
                %Meaning of the spectrums
                SPE=mean(SS(:,1:size(SS,2))');               
                %Calculation of the time vectors
                [Time_coh,~]=TIME(startTime,endTime,fs);
                %Concatenate the results
                SPE_COR(i)={SPE'};
                TIME_COH(i)={datenum(Time_coh)};
                FTOT(i)={F};
            elseif nny<params.n_sta
                %Error message
                set(text1, 'String', 'Error. Invaild Coherogram parameters.');
                drawnow
                return    
            else
                %Warning message
                set(text1, 'String', 'Warning. One or more stations are no available.');
                drawnow
            end
        end
        nny=size(SPE_COR,2);%Number of spectrums
        ZZ=[];%Initialize the maximum spectrum point vector
        FF=[];%Initialize the frequencies vector of the maximum spectral points
        %Loop through the number of spectrums
        for jj=1:1:nny
            %Find the maximum intensity of each spectrum and its frequency 
            [Cmax ncol]=max(SPE_COR{1,jj});
            Fmax=F(ncol);
            Z=max(SPE_COR{1,jj});
            %Concatenate the results
            FF=[FF;Fmax];
            ZZ=[ZZ;Z];
            Z=[];
        end
        %Set the position of the maximum spectral points 
        Z1=(ones(size(TIME_COH,1),1).*ZZ');
        Z1=Z1+Z1.*(5);
        %Plot the results in the coheregram axes
        ploting1(TIME_COH,SPE_COR,FTOT,FF,Z1,axes1)
        %Set the limits of x-axis on the coheregram and rms axes
        time=datetime(cell2mat(TIME_COH),'ConvertFrom','datenum');%Convert in datetime format
        if ~isempty(time) & ~isempty(cell2mat(TIME_AMP))
            xlim(axes1,[min(time) max(time)])
            xlim(axes3,[min(time) max(time)])
        elseif ~isempty(time) & isempty(cell2mat(TIME_AMP)) 
            xlim(axes1,[min(time) max(time)])
        end
        %Command message
        set(text1, 'String','Calculation finished.');
    catch
        %Error assessment
        if isempty(params.coord) | isempty(station);set(text1, 'String','Error. No coordinate file directory or No correct coordinate file.');drawnow;return;end
        if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
        if isempty(ff2);set(text1, 'String','Error. No files are found.');drawnow;return;end
        if params.f_high_pass>fs/2;set(text1, 'String','Error. Invalid High pass filter frequency.');drawnow;return;end
        if isempty(cell2mat(SPE_COR));set(text1, 'String', 'Error. Reset Coherogram parameters.');drawnow;return;end
        if isempty(z);set(text1, 'String','Error. No correct signals data.');drawnow;return;end
    end    
end
%==========================================================================
%%%%Calculation and plot of the rms%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calc2(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes3)%Clear the axes
    try
        %Acquisition and validation of the analysis parameters
        w_rms=str2num(get(h.w_rms,'string'));w_rms=round(w_rms,2);if w_rms<=0;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
        factor=str2num(get(h.factor,'string'));if factor<=0;set(text1, 'String','Error. Invalid temporal  factor mean.');drawnow;return;end   
        f1=str2num(get(h.f1,'string'));if f1<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        f2=str2num(get(h.f2,'string'));if f2<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        step_rms=str2num(get(h.step_rms,'string'));if step_rms<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        FREQ=f1:step_rms:f2;
        if isempty(FREQ) | length(FREQ)<2;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        %Command message
        textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
        set(text1, 'String', textLabel); 
        drawnow   
        %Reading .mat files
        [~, len]=size (ff);%Number of files selected
        TIME_AMP=cell(len,length(FREQ)-1);%Initialize the time matrix
        AMPP=cell(len,length(FREQ)-1);%Initialize the rms matrix
        %Loop through the frequency analysis
        for k=1:length(FREQ)-1
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
                amp=RMS(z,fs,dt,FREQ(k),FREQ(k+1));
                [~,Time_rms]=TIME(startTime,endTime,fs);
                %Concatenate the results
                TIME_RMS(i)={datenum(Time_rms)};
                AMP(i)={amp};
            end
            %Concatenate the results through the columns
            TIME_AMP(:,k)= TIME_RMS;
            AMPP(:,k)= AMP;
        end
        %Set the block of results
        slidval=length(FREQ)-1;
        %Plot the results in the rms axes
        ploting2(slidval,axes3)
        %Set the limits of x-axis on rms axes
        if ~isempty(cell2mat(TIME_COH))
            time=datetime(cell2mat(TIME_COH),'ConvertFrom','datenum');%Convert in datetime format
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
function  [time_coh,time_rms]=TIME(startTime,endTime,fs)
             time_coh=[];%Initialize the time vector of the coheregram
             time_rms=[];%Initialize the time vector of the rms
             %Convert the time limits of the trace in datetime format
             Time1=datetime(startTime,'ConvertFrom','datenum');
             Time2=datetime(endTime,'ConvertFrom','datenum');
             %Compute the time vector of the spectrums
             Time_coh=Time2-seconds(fs);
             time_coh=[time_coh;Time_coh];
             %Compute the time vector of the rms
             Time_rms_1=Time1+seconds(w_rms);
             Time_rms_2=Time2;
             Time_rms=Time_rms_1:seconds(w_rms):Time_rms_2;
             time_rms=[time_rms;Time_rms'];    
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
        i1=0; i2=0; k=0;
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
        if f1>=fs/2 | f2>fs/2 | step_rms>=fs/2;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    end
end

%==========================================================================
%%%%Plot of the coheregram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting1(time,coh,Ftot,FF,Z1,ax)
    try
        cla(ax)%Clear the coheregram axes
        time=datetime(cell2mat(time),'ConvertFrom','datenum');%Convert in datetime format
        coh=cell2mat(coh);
        Ftot=cell2mat(Ftot);
        %Logarithmic spectrums
        if logval==1
            coh=log10(coh);
            Z1=log10(Z1);
        end
        %Plot the coheregram
        surf(ax,time,Ftot,coh);
        shading(ax,'interp')
        axis(ax,'tight')
        hold(ax,'on')
        %Plot the maximum spectral points 
        plot3(ax,time,FF,Z1,'. k','MarkerSize',14);
        hold(ax,'off')
        ylim(ax,[0 15])
        ylabel(ax,'Frequency(Hz)');
        colormap(ax,'parula')
        colorbar(ax,'northoutside');
%         caxis(ax,[0 1])
        view ([0 90])
    catch
        %Error assessment
        if size(time,1)~=size(Ftot,1) & size(time,1)~=size(coh,1) & size(Ftot,1)~=size(coh,1)
            set(text1, 'String','Error. Incoherent dimensions of the Coherogram.');drawnow;return;
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
        plot(time,rms_mean./max(rms_mean))
        ylabel(ax,'Rms (a.u.)');
        legend(ax,[num2str(FREQ(ii)) '-' num2str(FREQ(ii+1)) ' Hz'])
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
   if isempty(cell2mat(SPE_COR));set(text1, 'String','Error. Coherogram matrix is empty.');drawnow;return;end
    if (get(h,'Value') == get(h,'Max'))
        %Command message
        set(text1, 'String', 'Logarithmic spectrum values.');
        %Set the logarithmic values of the spectrums
        logval=1;
        ploting1(TIME_COH,SPE_COR,FTOT,FF,Z1,axes1)
    else
        %Command message
        set(text1, 'String', 'Linear spectrum values.');
        %Set the logarithmic values of the spectrums
        logval=0;
        ploting1(TIME_COH,SPE_COR,FTOT,FF,Z1,axes1)
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
   slidval=slidval+1;if slidval>length(FREQ)-1;slidval=length(FREQ)-1;end
   %Plot the new block of rms
   ploting2(slidval,axes3);        
end

%==========================================================================
%%%%Set the limits of the y-axis on the coheregram axes%%%%%%%%%%%%%%%%%%%%
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
%%%%Set the limits of the x-axis on the coheregram and rms axes%%%%%%%%%%%%
function setlimits2(~,~)
    try
        %Manual settings
        %Acquisition time limits
        t1_lim=datetime(get(h.t1_lim,'string'),'InputFormat','yyMMdd-HHmmSS');%Convert in datetime format
        t2_lim=datetime(get(h.t2_lim,'string'),'InputFormat','yyMMdd-HHmmSS');%Convert in datetime format 
        xlim(axes1,[t1_lim t2_lim])
        xlim(axes3,[t1_lim t2_lim])
    catch
        %Default settings
        if ~isempty(cell2mat(TIME_COH)) 
            time=datetime(cell2mat(TIME_COH),'ConvertFrom','datenum');%Convert in datetime format
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
%%%%Save coheregram plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save1(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    try
        time=datetime(cell2mat(TIME_COH),'ConvertFrom','datenum');%Convert in datetime format
        coh=cell2mat(SPE_COR);
        Ftot=cell2mat(FTOT);
        %Error assessment
        if isempty(time) | isempty (coh);set(text1, 'String','Error. Coherogram matrix is empty.');drawnow;return;end
        %Plot the coheregram
        new_figure=figure('Name','Coherogram plot','NumberTitle','off');
        %Logarithmic spectrums
        if logval==1
            coh=log10(coh);
        end
        surf(time,Ftot,coh)
        shading('interp')
        axis('tight')
        hold('on')
        %Plot the maximum spectral points 
        plot3(time,FF,Z1,'. k','MarkerSize',14);
        hold('off')
        ylim([0 15])
        ylabel('Frequency(Hz)');
        colormap('parula')
        colorbar('eastoutside');
%         caxis([0 1])
        view ([0 90])
    catch
        %Error assessment
        if size(time,1)~=size(Ftot,1) & size(time,1)~=size(coh,1) & size(Ftot,1)~=size(coh,1)
            set(text1, 'String','Error. Incoherent dimensions of the Coherogram.');drawnow;return;
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
        legend([num2str(FREQ(ii)) '-' num2str(FREQ(ii+1)) ' Hz'])
        ylim([0 1])
    catch
        %Error assessment
        if isempty(time) | isempty (amp);set(text1, 'String','Error. Rms matrix is empty.');drawnow;return;end
        if size(amp,1)~=size(time,1);set(text1, 'String','Error. Incoherent dimensions of the Rms array.');drawnow;return;end
    end            
end

%==========================================================================
%%%%Save the results of the coheregram%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savespec(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Validation of the results
    if isempty(params.save);set(text1, 'String', sprintf('Error. No output directory.'));drawnow;return;end  
    if isempty(cell2mat(SPE_COR));set(text1, 'String', sprintf('Error. Coherogram matrix is empty.'));drawnow;return;end 
    if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
    %Loop through the number of Input files (ff)
    for kk=1:size(SPE_COR,2)
        time=TIME_COH{kk};spec=SPE_COR{kk};Time=time.*(ones(length(spec),1));
        Freq=FTOT{kk};name=ff{kk};[~,name,~] = fileparts(name);
        %Create the Output filename
        switch params.comp_chan 
            case 'Chan'
                name=name(end-14:end);
            case 'Comp'  
                name=name(end-18:end);
        end
        %Compute the juliand day
        jday=day(datetime(time(1),'ConvertFrom','datenum'),'dayofyear');
        if jday<10;jday=['00' num2str(jday)];elseif jday>=10 & jday<100;jday=['0' num2str(jday)];elseif jday>=100;jday=[num2str(jday)];end
        %Create the Output directory
        pathout=[params.save '/COH/' jday '/'];if ~exist(pathout, 'dir');mkdir(pathout);end 
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
    set(text1, 'String', sprintf(['COH' params.outputformat '-Observed data are saved.']));
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
            freq1=FREQ(jj);freq2=FREQ(jj+1);Ffolder=[num2str(freq1) '-' num2str(freq2) 'Hz'];
            name=ff{kk};name=name(end-27:end-4); 
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