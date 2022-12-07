%%%%%%%%%%%%% Calculation of the Polarization attributes %%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params: user parameters
% avfig: Polarization figure
% panel_color: colors of the main panel 
% fig_color: colors of the figure
% pan1: control right side panel 
% text1: Window message
% axes1,axes2,axes3: axes of the current figure
% az1_lim, az2_lim: limit of y-axis on azimuth axes
% inc1_lim, inc2_lim: limit of y-axis on incidence axes
% ret1_lim, ret2_lim: limit of y-axis on rectilinearity axes
% t1_lim, t2_lim: limit of x-axis on polarization axes
% slidval: value of dynamic plots 
% h: UI control parameters
%%%%%%%%%%%h.variables/params.variables/variables%%%%%%%%%%%%%%%%%%%%%%%%%%
% rad,Rad: radio button within a button group 
% rad1, rad2, rad3, rad4, Rad1, Rad2: button groups 
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
% list: list of 3-component seismic traces directories
% mytrace: seismic trace structure
% dt: sampling interval
% fs: sampling rate
% nyq: nyquist number
% z: amplitude vector of the seismic trace
% nt: number of samples
% startTime, endTime: temporal limits of the traces 
%--- Polarization setting
% typeMethod: type of analysis (Continous/Discrete)
% w_pol: the analysis window for the polarization analysis (s)
% f1_pol,f1: the minimum frequency for the polarization analysis (Hz)
% f2_pol,f1: the maximum frequency for the polarization analysis (Hz)
% step_pol: the frequency step for the polarization analysis (Hz)
% factor_pol,factor: the temporal factor for the moving average of the polarization attributes (hours)
% thr: the threshold value of rectlinearity 
% F: range of frequency analysis (Hz) 
% amp: amplitude vectors of the rms
% Time_pol, T_pol, TIME_POL: time vectors
% az, AZ,Azimuth: azimuth vectors/matrixs (degrees)
% incid,INC, Incidence: incidence vectors/matrixs (degrees)
% rect,RET, Rect: rectilinearity vectors/matrixs
% LL: sum of the auto-variances  
% TL: time samples


function polariz()
%% Declaration global parameters
clc
global params
comp=params.component;typeMethod='Con';
component=[];
ff= cell(1);
mytrace=[];
TIME_POL=[];
Azimuth=[];
Incidence=[];
Rect=[];
slidval=1;
pathname=[];
pathnam=[];
pathname2=[];
F=[];
f1=[];f2=[];step_pol=[];factor=[];w_pol=[];

%% Polarization figure
fig_color= [.94 .94 .94];
panel_color= [.97 .97 .97];

avfig= figure('Name', 'Polaritazion', 'position',get(0,'screensize'), ...
    'NumberTitle','off','units','normalized', 'color', fig_color,'MenuBar','none');

%Window message to user
text1= uicontrol('Style','edit','FontName','Arial','Fontsize',10,'Fontweight','bold',...
                'String','Window Message','Max',5,'Min',0,'HorizontalAlignment','left',...
                'BackgroundColor', 'w','Enable','inactive',...
                'units','normalized','Position',[.05 .02 .72 .05]);


%% Axeses
%Azimuth axes
axes1= axes('position',[.05 .72 .72 .25]);
hold(axes1,'on')
     grid(axes1,'on')
     grid(axes1,'minor')
     ylabel(axes1,'Azimuth (degrees)');
     set(axes1,'XTickLabel',[]);
     
%Incidence axes
axes2= axes('position',[.05 .42 .72 .25]);
hold(axes2,'on')
     grid(axes2,'on')
     grid(axes2,'minor')
     ylabel(axes2,'Incidence (degrees)');
     set(axes2,'XTickLabel',[]);

%Rectilinearity axes
axes3= axes('position',[.05 .12 .72 .25]);
hold(axes3,'on')
     grid(axes3,'on')
     grid(axes3,'minor')
     ylabel(axes3,'Rectilinearity');

%% Control right side panel
pan1= uipanel(avfig,'visible','on','Position',[.81 .03 .185 .94],...
    'BackgroundColor', panel_color);

uipanel(avfig,'visible','on','Position',[.81 .30 .185 .01],...
    'BackgroundColor', fig_color);


%% Text boxes
%-------------------------------- Polarization setting     
uicontrol(pan1,'Style','text', 'String','Pola setting',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.55 .71 .4 .06],...
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
uicontrol(pan1,'Style','text', 'String','Thr',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .43 .2 .04],...
    'BackgroundColor',panel_color)

%-------------------------------- Axis setting 
uicontrol(pan1,'Style','text', 'String','Axis setting',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.03 .73 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Az1 (deg)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .70 .4 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Az2 (deg)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .65 .4 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Inc1 (deg)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .58 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Inc2 (deg)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .53 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Ret1',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .48 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Ret2',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .43 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','t1',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .38 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','t2',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .33 .4 .04],...
    'BackgroundColor',panel_color);

%% Editable text 
%-------------------------------- Polarization setting
h.w_pol= uicontrol(pan1,'Style','edit', 'String',params.w_pol,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .69 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.f1= uicontrol(pan1,'Style','edit', 'String',params.f1_pol,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .64 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.f2= uicontrol(pan1,'Style','edit', 'String',params.f2_pol,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .59 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.step_pol= uicontrol(pan1,'Style','edit', 'String',params.step_pol,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .54 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.factor= uicontrol(pan1,'Style','edit', 'String',params.factor_pol,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .49 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.thr= uicontrol(pan1,'Style','edit', 'String','0.8',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .44 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);

%-------------------------------- Axis setting
h.az1_lim= uicontrol(pan1,'Style','edit', 'String','0',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.25 .69 .2 .04],...
    'BackgroundColor','w','callback',@setlimits1);
h.az2_lim= uicontrol(pan1,'Style','edit', 'String','360',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.25 .64 .2 .04],...
    'BackgroundColor','w','callback',@setlimits1);
h.inc1_lim= uicontrol(pan1,'Style','edit', 'String','0',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.25 .59 .2 .04],...
    'BackgroundColor','w','callback',@setlimits2);
h.inc2_lim= uicontrol(pan1,'Style','edit', 'String','90',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.25 .54 .2 .04],...
    'BackgroundColor','w','callback',@setlimits2);
h.ret1_lim= uicontrol(pan1,'Style','edit', 'String','0',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.25 .49 .2 .04],...
    'BackgroundColor','w','callback',@setlimits3);
h.ret2_lim= uicontrol(pan1,'Style','edit', 'String','1',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.25 .44 .2 .04],...
    'BackgroundColor','w','callback',@setlimits3);
h.t1_lim= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.25 .39 .3 .04],...
    'BackgroundColor','w','callback',@setlimits4);
h.t2_lim= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.25 .34 .3 .04],...
    'BackgroundColor','w','callback',@setlimits4);

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
uicontrol(pan1,'Style','text', 'String','Mode',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.63 .85 .4 .05],...
    'BackgroundColor',panel_color);

% Analysis Mode selection radio buttons
h.Rad= uibuttongroup(pan1,'units','normalized','BackgroundColor',panel_color,...
    'bordertype','none','Position',[.63 .82 .9 .05]);

set(h.Rad,'SelectionChangeFcn',@radcbk2);

h.Rad1 = uicontrol( h.Rad, 'Style','Radio','String','Con',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.0 0 .2 1],'HandleVisibility','off');
h.Rad2 = uicontrol( h.Rad, 'Style','Radio','String','Dis',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.21 0 .2 1],'HandleVisibility','off');

%Set the analysis mode (continous/discrete)
if strcmp(typeMethod,'Con')
    set (h.Rad1,'value',1);
    set (h.Rad2,'value',0);
elseif strcmp(typeMethod,'Dis')
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
                'Value',0,'Position',[.815 .32 .05 .03],'callback',@logscale);
% Editable text
uicontrol(pan1,'Style','text', 'String','Polarization',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.03 .21 .4 .07],...
    'BackgroundColor',panel_color);
% Calculation and save button
uicontrol('style','pushbutton', 'string','Calculation','units','normalized',...
    'position',[.815 .21 .08 .0494], 'callback',@calc); 
uicontrol('style','pushbutton', 'string','Save Plot','units','normalized',...
    'position',[.815 .13 .08 .0494], 'callback',@save1); 
uicontrol('style','pushbutton', 'string','Save Data','units','normalized',...
    'position',[.815 .05 .08 .0494], 'callback',@savepol);
uicontrol('style','pushbutton', 'string','Polar Scatter','units','normalized',...
    'position',[.905 .21 .08 .0494], 'callback',@Scatterpol);
uicontrol('style','pushbutton', 'string','Polar Hist','units','normalized',...
    'position',[.905 .13 .08 .0494], 'callback',@Histpol);

% Next and pre button
uicontrol('style','pushbutton', 'string','Next','units','normalized',...
    'position',[.776 .5 .03 .04], 'callback',@pola2);
uicontrol('style','pushbutton', 'string','Pre','units','normalized',...
    'position',[.776 .56 .03 .04], 'callback',@pola1);

%% Other functions (nested)
%%%%Return the list of files manually selected%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function openfnc(~,~)
    TIME_POL=[];Azimuth=[];Incidence=[];Rect=[];%Initialize Polarization attributes vectors
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);cla(axes2);cla(axes3)%Clear the axes
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
     TIME_POL=[];Azimuth=[];Incidence=[];Rect=[];%Initialize Polarization attributes vectors
     text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);cla(axes2);cla(axes3)%Clear the axes
    pathname=[];
    pathnam= uigetdir(params.data);
    pathnam=[pathnam '/'];
    search_sac=strcat(params.sta_ref,'.BH',comp,'*.mat');
    component=comp;
    
    if pathnam~=0
    % Get all files in the directory and subdirectories
    [ff,pathname2] = getAllFiles(pathnam,search_sac);
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
%%%%Return the type of analysis method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radcbk2(~,eventdata)
    
    typeMethod = get(eventdata.NewValue,'String');   
    
end

%==========================================================================
%%%%Calculation and plots the Polarization attributes%%%%%%%%%%%%%%%%%%%%%%
function calc(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);cla(axes2);cla(axes3)%Clear the axes
    try
        %Acquisition and validation of the analysis parameters
        w_pol=str2num(get(h.w_pol,'string'));w_pol=round(w_pol,2);if w_pol<=0;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
        factor=str2num(get(h.factor,'string'));if factor<=0;set(text1, 'String','Error. Invalid temporal  factor mean.');drawnow;return;end     
        f1=str2num(get(h.f1,'string'));if f1<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        f2=str2num(get(h.f2,'string'));if f2<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        step_pol=str2num(get(h.step_pol,'string'));if step_pol<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        F=f1:step_pol:f2;if isempty(F)|length(F)<2;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end   
        %Command message
        textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
        set(text1, 'String', textLabel); 
        drawnow
        %%Execute one of two groups of routines (Component/Channel)
        switch params.comp_chan
            case 'Comp'
                %Reading .mat files
                [~, len]=size (ff);%Number of files selected
                TIME_POL=cell(len,length(F)-1);%Initialize the time matrix
                Azimuth=cell(len,length(F)-1);%Initialize the azimuth matrix
                Incidence=cell(len,length(F)-1);%Initialize the incidence matrix
                Rect=cell(len,length(F)-1);%Initialize the rectilinearity matrix
                %Loop through the frequency analysis
                for k=1:length(F)-1
                    T_POL=cell(len,1);%Initialize the time vector
                    AZ=cell(len,1);%Initialize the azimuth vector
                    INC=cell(len,1);%Initialize the incidence vector
                    RET=cell(len,1);%Initialize the rectilinearity vector
                    %Loop through the number of files selected
                    for i=1:len
                        file_name=ff{i};%return the i-th file
                        file_name=file_name(end-18:end);%Create the term of research
                        %Control and set the Input folder
                        if ~isempty(pathnam);pathname=pathname2{i};pathname=[pathname '/'];end
                        %Search sethe 3 components of the seismic traces
                        file1=strcat(pathname,'/',params.sta_ref,'.*.',file_name);
                        list=dir(file1);
                        %Control the calculation through the number of components
                        if length(list)==3
                            sis=[];%Initialize the 3D traces
                            %Loop through the three components
                            for kk=1:length(list)    
                                filez=strcat(pathname,'/',list(kk).name);    
                                load(filez);%load the kk-th component 
                                z=mytrace.data;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate; nyq=0.5/dt; nt=mytrace.sampleCount;
                                startTime=mytrace.startTime;endTime=mytrace.endTime;
                                %Rescale and remove the best straight-fit line
                                z=z-mean(z);
                                z=detrend(z);
                                %Concatenate the seismic traces through the columns
                                sis(:,length(list)+1-kk)=z;
                            end
                            %Apply the Jurkevics algorithm for the calculation of the polarization attributes
                            [az,incid,rect,LL,TL] = filter_taper(sis,fs,F(k),F(k+1),w_pol,w_pol);
                            %Calculation of the time vectors
                            Time_pol=TIME(startTime,endTime,fs);
                            %Activation/Deactivation of the analysis type
                            if typeMethod=='Dis'
                                %Calculation of the rms vectors
                                amp=RMS(sis(:,1),fs,dt,f1,f2,w_pol,w_pol);
                                %Selection of discrete values
                                [val,imax]=max(amp);
                                Time_pol=Time_pol(imax);az=az(imax);incid=incid(imax);rect=rect(imax);    
                            end
                            %Concatenate the results
                            T_POL(i)={datenum(Time_pol)};
                            AZ(i)={az'};
                            INC(i)={incid'};
                            RET(i)={rect'};
                        else
                            %Warning message
                            textLabel = sprintf('Warning. One or more three components files are no available');
                            set(text1, 'String', textLabel); 
                            drawnow
                        end
                    end
                    %Concatenate the results through the columns
                    Azimuth(:,k)=AZ;
                    Incidence(:,k)=INC;
                    Rect(:,k)=RET;
                    TIME_POL(:,k)=T_POL;
                end
                %Set the block of results
                slidval=length(F)-1;
                %Plot the results in the azimuth axes
                ploting1(slidval,axes1)
                %Plot the results in the incidence axes
                ploting2(slidval,axes2)
                %Plot the results in the rectilinearity axes
                ploting3(slidval,axes3)
            case 'Chan'
                %Error message
                textLabel = sprintf('Error. No three components files are found.');
                set(text1, 'String', textLabel);
                drawnow
                return
        end
        %Command message
        set(text1, 'String','Calculation finished.');
        drawnow
    catch
        %Error assessment
        if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
        if isempty(list);set(text1, 'String','Error. No three components files are found');drawnow;return;end
        if isempty(z);set(text1, 'String','Error. No correct signals data.');drawnow;return;end
        if length(z)<=round(w_pol*fs);set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end 
        if f1>=fs/2 | f2>fs/2 | step_pol>=fs;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    end
end
%==========================================================================
%%%%Calculation of the time vector%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  time_pol=TIME(startTime,endTime,fs);
             time_pol=[];%Initialize the time vector
             %Convert the time limits of the trace in datetime format
             Time1=datetime(startTime,'ConvertFrom','datenum');
             Time2=datetime(endTime,'ConvertFrom','datenum');
             %Compute the time vector of the  polarization analysis
             w=w_pol;
             Time_pol_1=Time1+seconds(w_pol);
             Time_pol_2=Time2;
             Time_pol=Time_pol_1:seconds(w):Time_pol_2;
             time_pol=[time_pol;Time_pol'];  
end

%==========================================================================
%%%%Calculation of the rms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function amp=RMS(z,fs,dt,f1,f2,wls,wss)
    %Filter the seismic traces with a band-pass
    fs=(1/dt);fn=fs/2;
    W(1,1)=f1/fn; W(1,2)=f2/fn;
    [a b]=butter(2,W);
    z=filter(a,b,z);
    %Move the analysis window through the seismic traces
    N=length(z);wl=round(wls/dt);ws=round(wss/dt);
    i1=0; i2=0;k=0;
    while i2 < N-wl || i2== N-wl
        k=k+1;
        i1=(1+(k-1)*ws); i2=i1+wl-1;%Limits of the analysis window
        y1=z(i1:i2).*tukeywin(wl);%tapered cosine window
        %Apply the rms algorithm
        amp(k)=sqrt(sum((y1.*conj(y1))/wl));
    end
    amp=amp';
end
%==========================================================================
%%%%Plot of the Azimuth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting1(ii,ax)
    try
        thr=str2num(get(h.thr,'string'));%the threshold value of rectlinearity
        cla(ax)%Clear the azimuth axes
        hold(ax,'on')
        Az=cell2mat(Azimuth);
        Time=datetime(cell2mat(TIME_POL),'ConvertFrom','datenum');%Convert in datetime format
        idx2=cell2mat(Rect)>thr;%Selection of the values over the threshold
        Az=Az(:,ii);%Select the ii-th block of the azimuth
%         %Rescale the azimuth values under 180°
%         idx=Az<=180;
%         Az(idx)= Az(idx)+180;
        Time=Time(:,ii);%Select the ii-th block of the time
        idx2=idx2(:,ii);%Select the ii-th block of the selection
        Az=Az(idx2);%Azimuthal selection of the values over the threshold
        Time=Time(idx2);%Temporal selection of the values over the threshold
        %Moving average of the azimuth according the factor value
        meanangle_deg = rad2deg(atan2(movmean(sind(Az),hours(factor),'SamplePoints',Time),movmean(cosd(Az),hours(factor),'SamplePoints',Time)));
        az_mean = mod(meanangle_deg ,360);
        %Moving standard deviation of the azimuth according the factor value 
        meanangle_deg = rad2deg(atan2(movstd(sind(Az),hours(factor),'SamplePoints',Time),movstd(cosd(Az),hours(factor),'SamplePoints',Time)));
        az_std = mod(meanangle_deg ,360);
        %Plot the values of the azimuth
        plot(ax,Time,az_mean,'b')
        ylim(ax,[0 360])
        %Compute and plot the standard deviation of the azimuth
        time2=[Time', fliplr(Time')];
        curve1=az_mean+ az_std;
        curve2=az_mean- az_std;
        inBetween = [curve1', fliplr(curve2')];
        fill(ax,time2, inBetween, 'b','FaceAlpha',0.2,'EdgeColor','None');
        ylabel(ax,'Azimuth (degrees)');
        legend(ax,[num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
        set(ax,'XTickLabel',[]);
        hold (ax,'off')
    catch
        %Error assessment
        if size(Az,1)~=size(Time,1);set(text1, 'String','Error. Incoherent dimensions of the Polarization array.');drawnow;return;end
    end
end

%==========================================================================
%%%%Plot of the Incidence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting2(ii,ax)
    try
        thr=str2num(get(h.thr,'string'));%the threshold value of rectlinearity
        cla(ax)%Clear the incidence axes
        hold(ax,'on')
        Inc=cell2mat(Incidence);
        Time=datetime(cell2mat(TIME_POL),'ConvertFrom','datenum');%Convert in datetime format
        idx=cell2mat(Rect)>thr;%Selection of the values over the threshold
        Inc=Inc(:,ii);%Select the ii-th block of the incidence
        Time=Time(:,ii);%Select the ii-th block of the time
        idx=idx(:,ii);%Select the ii-th block of the selection
        Inc=Inc(idx);%Incidence selection of the values over the threshold
        Time=Time(idx);%Temporal selection of the values over the threshold
        inc_mean=movmean(Inc,hours(factor),'SamplePoints',Time);%Moving average of the incidence according the factor value 
        inc_std= movstd(Inc,hours(factor),'SamplePoints',Time);%Moving standard deviation of the incidence according the factor value 
        %Plot the values of the incidence
        plot(ax,Time,inc_mean,'k')
        ylim(ax,[0 90])
        %Compute and plot the standard deviation of the incidence
        time2=[Time', fliplr(Time')];
        curve1=inc_mean+inc_std;
        curve2=inc_mean-inc_std;
        inBetween = [curve1', fliplr(curve2')];
        fill(ax,time2, inBetween, 'k','FaceAlpha',0.2,'EdgeColor','None');
        ylabel(ax,'Incidence (degrees)');
        legend(ax,[num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
        set(ax,'XTickLabel',[]);
        hold (ax,'off')
    catch
        %Error assessment
        if size(Inc,1)~=size(Time,1);set(text1, 'String','Error. Incoherent dimensions of the Polarization array.');drawnow;return;end
    end
end

%==========================================================================
%%%%Plot of the Rectlinearity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting3(ii,ax)
    try
        thr=str2num(get(h.thr,'string'));%the threshold value of rectlinearity
        cla(ax)%Clear the rectlinearity axes
        hold(ax,'on')
        Ret=cell2mat(Rect);
        Time=datetime(cell2mat(TIME_POL),'ConvertFrom','datenum');%Convert in datetime format
        idx=cell2mat(Rect)>thr;%Selection of the values over the threshold
        Ret=Ret(:,ii);%Select the ii-th block of the rectlinearity
        Time=Time(:,ii);%Select the ii-th block of the time
        idx=idx(:,ii);%Select the ii-th block of the selection
        Ret=Ret(idx);%Rectlinearity selection of the values over the threshold
        Time=Time(idx);%Temporal selection of the values over the threshold
        ret_mean=movmean(Ret,hours(factor),'SamplePoints',Time);%Moving average of the rectlinearity according the factor value 
        ret_std=movstd(Ret,hours(factor),'SamplePoints',Time);%Moving standard devation of the rectlinearity according the factor value 
        %Plot the values of the rectlinearity
        plot(ax,Time,ret_mean,'g')
        ylim(ax,[0 1])
        %Compute and plot the standard deviation of the rectlinearity
        time2=[Time', fliplr(Time')];
        curve1=ret_mean+ret_std;
        curve2=ret_mean-ret_std;
        inBetween = [curve1', fliplr(curve2')];
        fill(ax,time2, inBetween, 'g','FaceAlpha',0.2,'EdgeColor','None');
        ylabel(ax,'Rectilinearity ');
        legend(ax,[num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
        hold (ax,'off')
    catch
        %Error assessment
        if size(Ret,1)~=size(Time,1);set(text1, 'String','Error. Incoherent dimensions of the Polarization array.');drawnow;return;end
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
        set(axes1,'yScale','log'); 
        set(axes2,'yScale','log');
        set(axes3,'yScale','log');
    else
        %Command message
        set(text1, 'String', 'Linear axis chosen.');
        %Set the linear scale
        set(axes1,'yScale','linear');
        set(axes2,'yScale','linear');
        set(axes3,'yScale','linear');
    end
end

%==========================================================================
%%%%Refresh the polarization axes through the pre button%%%%%%%%%%%%%%%%%%%
function pola1(~,~)
   %Error control
   if isempty(cell2mat(Azimuth));set(text1, 'String','Error. Polarization matrix is empty.');return;end
   %Decrease the slidval value
   slidval=slidval-1;if slidval<=0;slidval=1;end
   %Plot the new block of azimuth
   ploting1(slidval,axes1);
   %Plot the new block of incidence
   ploting2(slidval,axes2);
   %Plot the new block of rectilinearity
   ploting3(slidval,axes3);
        
end

%==========================================================================
%%%%Refresh the polarization axes through the next button%%%%%%%%%%%%%%%%%%
function pola2(~,~)
    %Error control
    if isempty(cell2mat(Azimuth));set(text1, 'String','Error. Polarization matrix is empty.');return;end
     %Increase the slidval value
    slidval=slidval+1;if slidval>length(F)-1;slidval=length(F)-1;end
    %Plot the new block of azimuth
    ploting1(slidval,axes1);
    %Plot the new block of incidence
    ploting2(slidval,axes2);
    %Plot the new block of rectilinearity
    ploting3(slidval,axes3);
        
end
%==========================================================================
%%%%Set the limits of the y-axis on the azimuth axes%%%%%%%%%%%%%%%%%%%%%%%
function setlimits1(~,~)
    try
        %Manual settings
        %Acquisition degrees limits
        az1_lim=str2num(get(h.az1_lim,'string'));
        az2_lim=str2num(get(h.az2_lim,'string')); 
        ylim(axes1,[az1_lim az2_lim])
    catch
        %Default settings
        ylim(axes1,[0 360])
        %Reset azimuth limits
        set(h.az1_lim,'string','0');
        set(h.az2_lim,'string','360');
    end
end

%==========================================================================
%%%%Set the limits of the y-axis on the incidence axes%%%%%%%%%%%%%%%%%%%%%
function setlimits2(~,~)
    try
        %Manual settings
        %Acquisition degrees limits
        inc1_lim=str2num(get(h.inc1_lim,'string'));
        inc2_lim=str2num(get(h.inc2_lim,'string')); 
        ylim(axes2,[inc1_lim inc2_lim])
    catch
        %Default settings
        ylim(axes2,[0 90])
        %Reset incidence limits
        set(h.inc1_lim,'string','0');
        set(h.inc2_lim,'string','90');
    end
end

%==========================================================================
%%%%Set the limits of the y-axis on the rectilinearity axes%%%%%%%%%%%%%%%%
function setlimits3(~,~)
    try
        %Manual settings
        %Acquisition limits
        ret1_lim=str2num(get(h.ret1_lim,'string'));
        ret2_lim=str2num(get(h.ret2_lim,'string'));  
        ylim(axes3,[ret1_lim ret2_lim])
    catch 
        %Default settings
        ylim(axes3,[0 1])
        %Reset rectilinearity limits
        set(h.ret1_lim,'string','0');
        set(h.ret2_lim,'string','1');
    end
end

%==========================================================================
%%%%Set the limits of the x-axis on the polarization axes%%%%%%%%%%%%%%%%%%
function setlimits4(~,~)
    try
        %Manual settings
        %Acquisition time limits
        t1_lim=datetime(get(h.t1_lim,'string'),'InputFormat','yyMMdd-HHmmSS');%Convert in datetime format
        t2_lim=datetime(get(h.t2_lim,'string'),'InputFormat','yyMMdd-HHmmSS');%Convert in datetime format 
        xlim(axes1,[t1_lim t2_lim])
        xlim(axes2,[t1_lim t2_lim])
        xlim(axes3,[t1_lim t2_lim])
    catch
        %Default settings
        if ~isempty(cell2mat(TIME_POL)) 
            ii=slidval;
            Time=datetime(cell2mat(TIME_POL),'ConvertFrom','datenum');%Convert in datetime format
            Time=Time(:,ii);%Select the ii-th block of the time
            xlim(axes1,[min(Time) max(Time)])
            xlim(axes2,[min(Time) max(Time)])
            xlim(axes3,[min(Time) max(Time)])
        end
    end
end

%==========================================================================
%%%%Control and reset automatically analysis parameters%%%%%%%%%%%%%%%%%%%%
function deflimits(~,~)
    %Acquisition analysis parameters
    w_pol=str2num(get(h.w_pol,'string'));
    factor=str2num(get(h.factor,'string'));   
    f1=str2num(get(h.f1,'string'));
    f2=str2num(get(h.f2,'string'));
    step_pol=str2num(get(h.step_pol,'string'));
    thr=str2num(get(h.thr,'string'));
    
    %Control and eventually reset parameters
    if (size(w_pol,1)==0);set(h.w_pol,'string',params.w_pol);end
    if (size(factor,1)==0);set(h.factor,'string',params.factor_pol);end
    if (size(f1,1)==0);set(h.f1,'string',params.f1_pol);end
    if (size(f2,1)==0);set(h.f2,'string',params.f2_pol);end
    if (size(step_pol,1)==0);set(h.step_pol,'string',params.step_pol);end
    if (size(thr,1)==0);set(h.thr,'string','0.8');end
end

%==========================================================================
%%%%Save polarization plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save1(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    thr=str2num(get(h.thr,'string'));%the threshold value of rectlinearity
    try
        ii=slidval;
        Az=cell2mat(Azimuth);
        Inc=cell2mat(Incidence);
        ret=cell2mat(Rect);
        Time=datetime(cell2mat(TIME_POL),'ConvertFrom','datenum');%Convert in datetime format
        idx2=cell2mat(Rect)>thr;%Selection of the values over the threshold
        Az=Az(:,ii);%Select the ii-th block of the azimuth
%         %Rescale the azimuth values under 180°
%         idx=Az<=180;
%         Az(idx)= Az(idx)+180;
        Inc=Inc(:,ii);%Select the ii-th block of the incidence
        ret=ret(:,ii);%Select the ii-th block of the rectilinearity
        Time=Time(:,ii);%Select the ii-th block of the time
        idx2=idx2(:,ii);%Select the ii-th block of the selection 
        Az=Az(idx2);%Azimuthal selection of the values over the threshold
        Inc=Inc(idx2);%Incidence selection of the values over the threshold
        ret=ret(idx2);%Rectilinearity selection of the values over the threshold
        Time=Time(idx2);%temporal selection of the values over the threshold
        %Moving average of the azimuth according the factor value
        meanangle_deg = rad2deg(atan2(movmean(sind(Az),hours(factor),'SamplePoints',Time),movmean(cosd(Az),hours(factor),'SamplePoints',Time)));
        az_mean = mod(meanangle_deg ,360);
        %Moving standard deviation of the azimuth according the factor value 
        meanangle_deg = rad2deg(atan2(movstd(sind(Az),hours(factor),'SamplePoints',Time),movstd(cosd(Az),hours(factor),'SamplePoints',Time)));
        az_std = mod(meanangle_deg ,360);
        inc_mean=movmean(Inc,hours(factor),'SamplePoints',Time);%Moving standard devation of the incidence according the factor value 
        ret_mean=movmean(ret,hours(factor),'SamplePoints',Time);%Moving standard devation of the rectlinearity according the factor value   
        inc_std=movstd(Inc,hours(factor),'SamplePoints',Time);%Moving standard devation of the incidence according the factor value 
        ret_std=movstd(ret,hours(factor),'SamplePoints',Time);%Moving standard devation of the rectlinearity according the factor value 
        %Compute the standard deviation of the azimuth
        time2=[Time', fliplr(Time')];
        curveaz1=az_mean+az_std;
        curveaz2=az_mean-az_std;
        inBetweenaz = [curveaz1', fliplr(curveaz2')];
         %Compute the standard deviation of the incidence
        curveinc1=inc_mean+inc_std;
        curveinc2=inc_mean-inc_std;
        inBetweeninc = [curveinc1', fliplr(curveinc2')];
        %Compute the standard deviation of the rectlinearity
        curveret1=ret_mean+ret_std;
        curveret2=ret_mean-ret_std;
        inBetweenret = [curveret1', fliplr(curveret2')];
        new_figure=figure('Name','Polarization Attributes plot','NumberTitle','off');
        %Plot the values of the azimuth
        subplot(3,1,1)
        plot(Time,az_mean,'b')
        ylim([0 360])
        %Plot the standard deviation of the azimuth
        hold on
        fill(time2, inBetweenaz, 'b','FaceAlpha',0.2,'EdgeColor','None');
        hold off
        ylabel('Azimuth (degrees)');
        grid on
        grid minor
        legend([num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
        %Plot the values of the incidence
        subplot(3,1,2)
        plot(Time,inc_mean,'k')
        ylim([0 90])
        %Plot the standard deviation of the incidence
        hold on
        fill(time2, inBetweeninc, 'k','FaceAlpha',0.2,'EdgeColor','None');
        hold off
        ylabel('Incidence (degrees)');
        grid on
        grid minor
        legend([num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
        %Plot the values of the rectlinearity
        subplot(3,1,3)
        plot(Time,ret_mean,'g')
        ylim([0 1])
        %Plot the standard deviation of the rectilinearity
        hold on
        fill(time2, inBetweenret, 'g','FaceAlpha',0.2,'EdgeColor','None');
        hold off
        ylabel('Rectilinearity');
        grid on
        grid minor
        legend([num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
    catch
        %Error assessment
        if isempty(cell2mat(Time)) | isempty (cell2mat(Az));set(text1, 'String','Error. Polarization matrix is empty.');drawnow;return;end
        if size(Az,1)~=size(Time,1);set(text1, 'String','Error. Incoherent dimensions of the Polarization array.');drawnow;return;end
    end            
end

%==========================================================================
%%%%%%%%Polar scatter plot of Polarization Attributes%%%%%%%%%%%%%%%%%%%%%%
function Scatterpol(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    thr=str2num(get(h.thr,'string'));%the threshold value of rectlinearity
    try
        ii=slidval;
        Az=cell2mat(Azimuth);
        Inc=cell2mat(Incidence);
        ret=cell2mat(Rect);
        Time=datetime(cell2mat(TIME_POL),'ConvertFrom','datenum');%Convert in datetime format
        idx2=cell2mat(Rect)>thr;%Selection of the values over the threshold
        Az=Az(:,ii);%Select the ii-th block of the azimuth
        Inc=Inc(:,ii);%Select the ii-th block of the incidence
        ret=ret(:,ii);%Select the ii-th block of the rectilinearity
        Time=Time(:,ii);%Select the ii-th block of the time
        idx2=idx2(:,ii);%Select the ii-th block of the selection 
        Az=Az(idx2);%Azimuthal selection of the values over the threshold
        Inc=Inc(idx2);%Incidence selection of the values over the threshold
        ret=ret(idx2);%Rectilinearity selection of the values over the threshold
        Time=Time(idx2);%temporal selection of the values over the threshold       
        new_figure=figure('Name','Polar Scatter of Polarization Attributes','NumberTitle','off');
        %Expressing azimuth in radians
        ros=Az;
        fin=(2*pi*ros)./360;
        bins2=[0:5*((2*pi)/360):355*((2*pi)/360)];
        %Plotting polar scatter of azimuth, incidence and rectilinearity
        polarscatter(fin,Inc,[],ret,'filled','MarkerFaceAlpha',.2);
        set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
        title(['Azimuth (degrees)' '-' num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
        thetaticklabels({'N','30°','60°','E','120°','150°','S','210°','240°','W','300°','330°'})
        colorbar
        caxis([0 1])
        rlim([0 90])
    catch 
        %Error assessment
        if isempty(Time) | isempty(Az);set(text1, 'String','Error. Polarization matrix is empty.');drawnow;return;end
        if size(Az,1)~=size(Time,1);set(text1, 'String','Error. Incoherent dimensions of the Polarization array.');drawnow;return;end
    end
end

%==========================================================================
%%%%%%%%Polar histogram plot of Azimuth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Histpol(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    thr=str2num(get(h.thr,'string'));%the threshold value of rectlinearity
    try
        ii=slidval;
        Az=cell2mat(Azimuth);
        idx2=cell2mat(Rect)>thr;%Selection of the values over the threshold
        Az=Az(:,ii);%Select the ii-th block of the azimuth
        idx2=idx2(:,ii);%Select the ii-th block of the selection 
        Az=Az(idx2);%Azimuthal selection of the values over the threshold
        new_figure=figure('Name','Polar histogram of Azimuth','NumberTitle','off');
        %Expressing azimuth in radians
        ros=Az;
        fin=(2*pi*ros)./360;
        bins2=[0:5*((2*pi)/360):355*((2*pi)/360)];
        %Plotting polar histogram of azimuth
        newh=polarhistogram(fin,bins2,'FaceColor','b');
        newh.Normalization = 'probability';
        set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
        title(['Azimuth (degrees)' '-' num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
        thetaticklabels({'N','30°','60°','E','120°','150°','S','210°','240°','W','300°','330°'})
    catch 
        %Error assessment
        if isempty(Az);set(text1, 'String','Error. Polarization matrix is empty.');drawnow;return;end       
    end
end
%==========================================================================
%%%%Save the results of the polarization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savepol(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Validation of the results
    if isempty(params.save);set(text1, 'String', sprintf('Error. No output directory.'));drawnow;return;end  
    if isempty(cell2mat(Azimuth));set(text1, 'String', sprintf('Error. Polarization matrix is empty.'));drawnow;return;end
    if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
    %Loop through the frequency analysis
    for jj=1:size(Azimuth,2)
        %Loop through the number of Input files (ff)
        for kk=1:size(Azimuth,1)
            time=TIME_POL{kk,jj};az=Azimuth{kk,jj};inc=Incidence{kk,jj};ret=Rect{kk,jj}; 
            %Compute the juliand day
            jday=day(datetime(time(1),'ConvertFrom','datenum'),'dayofyear');
            if jday<10;jday=['00' num2str(jday)];elseif jday>=10 & jday<100;jday=['0' num2str(jday)];elseif jday>=100;jday=[num2str(jday)];end
            %Create the Output directory and filename
            freq1=F(jj);freq2=F(jj+1);Ffolder=[num2str(freq1) '-' num2str(freq2) 'Hz'];
            name=ff{kk};name=[params.sta_ref '.' name(end-18:end-4)];
            pathout=[params.save '/POL/' Ffolder '/' jday '/'];if ~exist(pathout, 'dir');mkdir(pathout);end 
            output=[pathout name params.outputformat];
            %Create a table of the results
            data=table(time,az,inc,ret);data.Properties.VariableNames={'Time','Azimuth','Incidence','Rectilinearity'};
            data.Properties.VariableUnits = {'datenum','degrees','degrees',' '};
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
    set(text1, 'String', sprintf(['POL' params.outputformat '-Observed data are saved.'])); 
    drawnow
end
end