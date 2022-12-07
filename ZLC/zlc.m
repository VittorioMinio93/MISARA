%%%%%%%%%%%%% Calculation of the ZLC parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params: user parameters
% avfig: ZLC figure
% panel_color: colors of the main panel 
% fig_color: colors of the figure
% pan1: control right side panel 
% text1: Window message
% axes1,axes2,axes3: axes of the current figure
% az1_lim, az2_lim: limit of y-axis on azimuth axes
% ray1_lim, ray2_lim: limit of y-axis on ray parameter axes
% cor1_lim, cor2_lim: limit of y-axis on cross correlation axes
% t1_lim, t2_lim: limit of x-axis on ZLC axes
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
% xs, ys, zs, utm_zone: coordinates in UTM system
% LAT, LON, ELE: coordinates in wgs84 system
%--- Data file
% data: the Input folder
% save: the Output folder
% outputformat: the format (.mat/.txt) of the output files
% search_sac, date_search: term of research file 
% pathname/pathnam/pathaname2, fname: Input folder and list of filenames
% ff: list of files with their path
% list: list of seismic traces directories
% mytrace: seismic trace structure
% dt: sampling interval
% fs: sampling rate
% nyq: nyquist number
% z: amplitude vector of the seismic trace
% nt: number of samples
% startTime, endTime: temporal limits of the traces 
%--- ZLC setting
% typeMethod: type of analysis (Continous/Discrete)
% emod: Error analysis routine (on/off)
% w_zlc: the analysis window for the ZLC analysis (s)
% f1_zlc,f1: the minimum frequency for the ZLC analysis (Hz)
% f2_zlc,f2: the maximum frequency for the ZLC analysis (Hz)
% step_zlc: the frequency step for the ZLC analysis (Hz)
% vel_zlc, v: the velocity value (km/s)
% splint: the Spline interpolation value (1 o 0)
% maxlags: the max lags for the Correlation coefficient calculation (s)
% F: range of frequency analysis (Hz)  
% amp: amplitude vectors of the rms
% thr: the threshold value of cross correlation 
% az1,az2,azstep: azimuthal limits and step for histograms 3D (degrees) 
% ray1,ray2,raystep: slowness limits and step for histograms 3D (s/km) 
% xcor1,xcor2,xcorstep: correlation limits and step for histograms 3D 
% t1,t2,nbins/timestep: temporal limits and step for histograms 3D
% az: azimuth vector
% baz,AZ,Azimuth: backazimuth vector/matrix  (degrees)
% rayp,RAYP,RayP: ray parameter vector/matrix (s/km)  
% delay: delay time computed by ZLC analysis (samples)
% xcmax,Xcorr,XCORR: cross correlation vector/matrix 
% inc,INC,Incidence: incidence vector/matrix (degrees)
% sx/SX/ssx, sy/SY/ssy: horizontal slowness vectors/matrixs (s/km) 
% ERROR1,Eb: backazimuth error vector/matrix (degrees)
% ERROR2,Er: ray parameter error vector/matrix (s/km)  
% ERROR3,Ei: incidence error vector/matrix (degrees)
% J_error: JackKnife results
% Time_zlc,T_zlc,TIME_ZLC: time vectors

function zlc()
%% Declaration global parameters
clc
global params
comp=params.component;typeMethod='Con';component=[];
ff= cell(1);
mytrace=[];
slidval=1;
pathnam=[];pathname=[];pathname2=[];
F=[];
Azimuth=[];
RayP=[];
TIME_ZLC=[];
XCORR=[];
Incidence=[];Eb=[];Er=[];Ei=[];
ssx=[];
ssy=[];
xs=[];ys=[];zs=[];station={};
LAT=[];LON=[];ELE=[];utm_zone=[];emod='n';

%% ZLC figure
fig_color= [.94 .94 .94];
panel_color= [.97 .97 .97];

avfig= figure('Name', 'Zero Lag CrossCorrelation', 'position',get(0,'screensize'), ...
    'NumberTitle','off','units','normalized', 'color', fig_color,'MenuBar','none');

%Window message to user
text1= uicontrol('Style','edit','FontName','Arial','Fontsize',10,'Fontweight','bold',...
                'String','Window Message','Max',5,'Min',0,'HorizontalAlignment','left',...
                'BackgroundColor', 'w','Enable','inactive',...
                'units','normalized','Position',[.05 .02 .72 .05]);


%% Axeses
% Cross Correlation axes
axes1= axes('position',[.05 .72 .72 .25]);
hold(axes1,'on')
     grid(axes1,'on')
     grid(axes1,'minor')
     ylim(axes1,[0 1])
     ylabel(axes1,'Correlation');
     colorbar(axes1)
     colormap(axes1,'parula')
% BackAzimuth axes
axes2= axes('position',[.05 .42 .72 .25]);
hold(axes2,'on')
     grid(axes2,'on')
     grid(axes2,'minor')
     ylim(axes2,[0 360])
     ylabel(axes2,'Back Azimuth (degrees)');
     colorbar(axes2)
     colormap(axes2,'parula')
% Rayparameter axes
axes3= axes('position',[.05 .12 .72 .25]);
hold(axes3,'on')
     grid(axes3,'on')
     grid(axes3,'minor')
     ylim(axes3,[0 1])
     ylabel(axes3,'RayParameter (s/km)');
     colorbar(axes3)
     colormap(axes3,'parula')

%% Control right side panel
pan1= uipanel(avfig,'visible','on','Position',[.81 .03 .185 .96],...
    'BackgroundColor', panel_color);

uipanel(avfig,'visible','on','Position',[.81 .30 .185 .01],...
    'BackgroundColor', fig_color);

%% Text boxes
%-------------------------------- ZLC setting 
uicontrol(pan1,'Style','text', 'String','ZLC setting',...
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
uicontrol(pan1,'Style','text', 'String','Vel (km/s)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .48 .3 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Xbin (h)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .43 .2 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Thr',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .38 .2 .04],...
    'BackgroundColor',panel_color);

%-------------------------------- Axis setting
uicontrol(pan1,'Style','text', 'String','Axis setting',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.03 .73 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Cor1',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .70 .4 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Cor2',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .65 .4 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Az1 (deg)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .58 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Az2 (deg)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .53 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Ray1 (s/km)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .48 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Ray2 (s/km)',...
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
%-------------------------------- ZLC setting
h.w_zlc= uicontrol(pan1,'Style','edit', 'String',params.w_zlc,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .69 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.f1= uicontrol(pan1,'Style','edit', 'String',params.f1_zlc,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .64 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.f2= uicontrol(pan1,'Style','edit', 'String',params.f2_zlc,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .59 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.step_zlc= uicontrol(pan1,'Style','edit', 'String',params.step_zlc,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .54 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.vel_zlc= uicontrol(pan1,'Style','edit', 'String',params.vel_zlc,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .49 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.nbins= uicontrol(pan1,'Style','edit', 'String','3',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .44 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.thr= uicontrol(pan1,'Style','edit', 'String','0.75',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .39 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);

%-------------------------------- Axis setting
h.cor1_lim= uicontrol(pan1,'Style','edit', 'String','0',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .69 .2 .04],...
    'BackgroundColor','w','callback',@setlimits1);
h.cor2_lim= uicontrol(pan1,'Style','edit', 'String','1',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .64 .2 .04],...
    'BackgroundColor','w','callback',@setlimits1);
h.az1_lim= uicontrol(pan1,'Style','edit', 'String','0',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .59 .2 .04],...
    'BackgroundColor','w','callback',@setlimits2);
h.az2_lim= uicontrol(pan1,'Style','edit', 'String','360',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .54 .2 .04],...
    'BackgroundColor','w','callback',@setlimits2);
h.ray1_lim= uicontrol(pan1,'Style','edit', 'String','0',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .49 .2 .04],...
    'BackgroundColor','w','callback',@setlimits3);
h.ray2_lim= uicontrol(pan1,'Style','edit', 'String','5',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .44 .2 .04],...
    'BackgroundColor','w','callback',@setlimits3);
h.t1_lim= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.183 .39 .3 .04],...
    'BackgroundColor','w','callback',@setlimits4);
h.t2_lim= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.183 .34 .3 .04],...
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
% Error mode checkbox           
uicontrol('Style','checkbox','BackgroundColor', panel_color,...
                'String','Error Mode','units','normalized','Tag','errmod',...
                'Value',0,'Position',[.925 .775 .06 .05],'callback',@errormode);
% Editable text
uicontrol(pan1,'Style','text', 'String','ZLC',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.03 .20 .4 .07],...
    'BackgroundColor',panel_color);
% Calculation and save button
uicontrol('style','pushbutton', 'string','Calculation','units','normalized',...
    'position',[.815 .2 .08 .0494], 'callback',@calc); 
uicontrol('style','pushbutton', 'string','Save Main Plot','units','normalized',...
    'position',[.815 .125 .08 .0494], 'callback',@save1); 
uicontrol('style','pushbutton', 'string','Save Data','units','normalized',...
    'position',[.815 .05 .08 .0494], 'callback',@save_zlc);
uicontrol('style','pushbutton', 'string','Ray vs. Az','units','normalized',...
    'position',[.905 .2 .08 .0494], 'callback',@save2);
uicontrol('style','pushbutton', 'string','Slow Grid','units','normalized',...
    'position',[.905 .125 .08 .0494], 'callback',@save3);
uicontrol('style','pushbutton', 'string','Hist','units','normalized',...
    'position',[.905 .05 .08 .0494], 'callback',@save4);


% Next and pre button
uicontrol('style','pushbutton', 'string','Next','units','normalized',...
    'position',[.776 .5 .03 .04], 'callback',@loca2);
uicontrol('style','pushbutton', 'string','Pre','units','normalized',...
    'position',[.776 .56 .03 .04], 'callback',@loca1);


%% Other functions (nested)
%%%%Return the list of files manually selected%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function openfnc(~,~)
     Azimuth=[];RayP=[];TIME_ZLC=[];XCORR=[];Incidence=[];Eb=[];Er=[];Ei=[];ssx=[];ssy=[];%Initialize ZLC attributes vectors
     xs=[];ys=[];zs=[];station={};LAT=[];LON=[];ELE=[];%Initialize coordinates vectors
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
    Azimuth=[];RayP=[];TIME_ZLC=[];XCORR=[];Incidence=[];Eb=[];Er=[];Ei=[];ssx=[];ssy=[];%Initialize ZLC attributes vectors
    xs=[];ys=[];zs=[];station={};LAT=[];LON=[];ELE=[];%Initialize coordinates vectors
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);cla(axes2);cla(axes3)%Clear the axes
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
%%%%Return the type of analysis method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radcbk2(~,eventdata)
    
    typeMethod = get(eventdata.NewValue,'String');   
    
end

%==========================================================================
%%%%Calculation and plots the ZLC attributes%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calc(~,~) 
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);cla(axes2);cla(axes3)%Clear the axes
    try
        %Acquisition and validation of the analysis parameters
        nbins=str2num(get(h.nbins,'string'));if nbins<=0;set(text1, 'String','Error. Invalid temporal bin value.');drawnow;return;end 
        v=str2num(get(h.vel_zlc,'string'));if v<=0;set(text1, 'String','Error. Invalid velocity value.');drawnow;return;end 
        splint=params.splint;if (splint~=0) & (splint~=1) ;set(text1, 'String','Error. Invalid Spline interpolation value.');drawnow;return;end
        maxlags=params.maxlags;if maxlags<=0;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
        w_zlc=str2num(get(h.w_zlc,'string'));if w_zlc<=0;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end  
        f1=str2num(get(h.f1,'string'));if f1<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        f2=str2num(get(h.f2,'string'));if f2<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        step_zlc=str2num(get(h.step_zlc,'string'));if step_zlc<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        F=f1:step_zlc:f2;if isempty(F)|length(F)<2;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        xs=[];ys=[];zs=[];station={};LAT=[];LON=[];ELE=[];%Initialize coordinates vectors
        %Execute one of two groups of routines (coordinates in wgs84/UTM)
        switch params.coord_sys
            case 'degrees'
                %Load station coordinate file
                load(params.coord);
                %Error control
                if isempty(LAT);set(text1, 'String','Error. No correct coordinate file.');drawnow;return;end
                %Convert WGS84 coordinates (Latitude, Longitude) into UTM coordinates
                [xs,ys,~]=wgs2utm(LAT,LON);
                zs=ELE;
                %Express coordinates in km 
                xs=xs./1000;
                ys=ys./1000;
                zs=zs./1000;
            case 'm'
                %Load station coordinate file
                load(params.coord)
                %Error control
                if isempty(xs);set(text1, 'String','Error. No correct coordinate file.');drawnow;return;end
        end
        X=[xs ys];%Coordinates matrix
        nstation=length(station);%Number of stations used
        %Command message
        textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
        set(text1, 'String', textLabel); 
        drawnow
        %Reading .mat files
        [~, len]=size (ff);%Number of files selected
        Azimuth=cell(len,length(F)-1);%Initialize the backazimuth matrix
        RayP=cell(len,length(F)-1);%Initialize the ray parameter matrix
        TIME_ZLC=cell(len,length(F)-1);%Initialize the time matrix
        XCORR=cell(len,length(F)-1);%Initialize the cross-correlation matrix
        ssx=cell(len,length(F)-1);ssy=cell(len,length(F)-1);%Initialize the horizontal slowness matrix
        Eb=cell(len,length(F)-1);Er=cell(len,length(F)-1);Ei=cell(len,length(F)-1);%Initialize the error matrixs
        Incidence=cell(len,length(F)-1);%Initialize the incidence matrix
        %Loop through the frequency analysis
        for kf=1:length(F)-1
            AZ=cell(len,1);%Initialize the backazimuth matrix
            RAYP=cell(len,1);%Initialize the ray parameter matrix
            SX=cell(len,1);SY=cell(len,1);%Initialize the horizontal slowness matrix
            T_zlc=cell(len,1);%Initialize the time matrix
            Xcorr=cell(len,1);%Initialize the cross-correlation matrix
            ERROR1=cell(len,1);ERROR2=cell(len,1);ERROR3=cell(len,1);%Initialize the error matrixs
            INC=cell(len,1);%Initialize the incidence matrix
            %Loop through the number of files selected
            for i=1:len
                file_name=ff{i};%return the i-th file
                file_name=file_name(end-18:end-4);%Create the term of research
                %Control and set the Input folder
                if ~isempty(pathnam);pathname=pathname2{i};pathname=[pathname '/'];end
                %Search all seismic traces based on the station names
                list=searchAllStaz(pathname,station,component,file_name,params.comp_chan);
                Z=[];%Initialize the signals matrix
                %Control the calculation through the number of stations
                if length(list)==nstation
                    %Loop through the number of stations used
                    for kk=1:length(list)    
                        filez=strcat(pathname,list(kk).name);    
                        load(filez);%load the ifile-th seismic trace  
                        z=mytrace.data;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate; nyq=0.5/dt; nt=mytrace.sampleCount;
                        startTime=mytrace.startTime;endTime=mytrace.endTime;
                        %Rescale and remove the best straight-fit line
                        z=z-mean(z);
                        z=detrend(z);
                        %Filter the seismic traces with a band-pass
                        fs=(1/dt);fn=fs/2;
                        f1=F(kf);
                        f2=F(kf+1);
                        W(1,1)=f1/fn; W(1,2)=f2/fn;
                        [a b]=butter(2,W);
                        z=filter(a,b,z);
                        %Concatenate the seismic traces through the columns
                        Z=[Z z];
                    end
                    N=length(Z);wl=round(w_zlc/dt);
                    ip=0;
                    %Loop through all independent pairs of stations
                    for is1=1:nstation-1
                        for is2=is1+1:nstation
                            ip=ip+1;
                            i1=0; i2=0;k=0;
                            %Move the analysis window through the seismic traces
                            while i2 < N-wl || i2== N-wl
                                k=k+1;
                                i1=(1+(k-1)*wl); i2=i1+wl-1;%Limits of the analysis window
                                y1=Z(i1:i2,is1).*tukeywin(wl); y2=Z(i1:i2,is2).*tukeywin(wl);%tapered cosine windows
                                %Calculation of the Cross-correlation coefficients and the delay times 
                                [delay(ip,k) xcmax(ip,k)]=GetDelay(y1,y2,splint,dt,maxlags);
                                %Apply the rms algorithm
                                if is1==1
                                    amp(k)=sqrt(sum((y1.*conj(y1))/wl));
                                end
                            end
                        end
                    end
                    %Calculation of the kinematic parameters and time vector
                    [sx sy az rayp]=FitPlane(X,-delay);
                    Time_zlc=TIME(startTime,endTime,fs);
                    %Execute one of two groups of routines (continous/discrete analysis)
                    switch typeMethod
                        case 'Dis'
                            %Selection of discrete values
                            xcmax=mean(xcmax)';[val,imax]=max(xcmax);xcmax=xcmax(imax);
                            Time_zlc=Time_zlc(imax);sx=sx(imax);sy=sy(imax);rayp=rayp(imax);az=az(imax);
                            %Convert azimuth in backazimuth
                            baz=az+180;
                            idx=baz>360;
                            %Activation/Deactivation of the error analysis
                            if emod=='y'
                                %Estimation of error analysis through JackKnife method
                                J_error=JackKnife_low(Z,splint,dt,maxlags,X,wl,wl,baz,rayp,sx,sy);
                                %Estimation of incidence errors
                                angl=(J_error(:,2).*v);
                                L=find(angl>1);
                                angl(L)=NaN;
                                %Concatenate the errors
                                ERROR1(i)={J_error(:,1)};
                                ERROR2(i)={J_error(:,2)};
                                ERROR3(i)={asind(angl)};
                            end
                            baz(idx)=baz(idx)-360;
                        case 'Con' 
                            %Convert azimuth in backazimuth
                            baz=az+180;
                            idx=baz>360;
                            %Activation/Deactivation of the error analysis
                            if emod=='y'
                                %Estimation of error analysis through JackKnife method
                                J_error=JackKnife_con(Z,splint,dt,maxlags,X,wl,wl,baz,rayp,sx,sy);
                                %Estimation of incidence errors
                                angl=(J_error(:,2).*v);
                                L=find(angl>1);
                                angl(L)=NaN;
                                %Concatenate the errors
                                ERROR1(i)={J_error(:,1)};
                                ERROR2(i)={J_error(:,2)};
                                ERROR3(i)={asind(angl)};
                            end
                            baz(idx)=baz(idx)-360;
                    end
                    %Estimation of incidence
                    angl=(rayp.*v);
                    L=find(angl>1);
                    angl(L)=NaN;
                    %Concatenate the results
                    INC(i)={asind(angl)};
                    AZ(i)={baz};
                    RAYP(i)={rayp};
                    SX(i)={sx};
                    SY(i)={sy};
                    T_zlc(i)={datenum(Time_zlc)};
                    Xcorr(i)={mean(xcmax)'};      
                else
                    %Warning message
                    set(text1, 'String', 'Warning. One or more stations are no available.');
                    drawnow
                end
            end
            %Error control
            if isempty(cell2mat(AZ));set(text1, 'String', 'Error. Reset Coordinates parameters.');drawnow;return;end
            %Concatenate the results through the columns
            Azimuth(:,kf)=AZ;
            RayP(:,kf)=RAYP;
            TIME_ZLC(:,kf)= T_zlc;
            XCORR(:,kf)= Xcorr;
            ssx(:,kf)=SX;
            ssy(:,kf)=SY;
            Eb(:,kf)=ERROR1;
            Er(:,kf)=ERROR2;
            Ei(:,kf)=ERROR3;
            Incidence(:,kf)=INC;
        end
        %Set the block of results
        slidval=length(F)-1;
        %Plot the results in the cross correlation axes
        ploting1(slidval,axes1)
        %Plot the results in the backazimuth axes
        ploting2(slidval,axes2)
        %Plot the results in the ray parameter axes
        ploting3(slidval,axes3)
        %Command message
        set(text1, 'String','Calculation finished.');
        drawnow
    catch
        %Error assessment
        if isempty(params.coord) | isempty(station);set(text1, 'String','Error. No coordinate file directory or No correct coordinate file.');drawnow;return;end
        if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
        if isempty(list);set(text1, 'String','Error. No files are found.');drawnow;return;end
        if isempty(z);set(text1, 'String','Error. No correct signals data.');drawnow;return;end
        if length(z)<=round(w_zlc*fs);set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
        if f1>=fs/2 | f2>fs/2 | step_zlc>=fs/2;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    end
end

%==========================================================================
%%%%Calculation of the time vector%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  time_zlc=TIME(startTime,endTime,fs)
             w_zlc=str2num(get(h.w_zlc,'string'));%Get the analysis window (s)
             time_zlc=[];%Initialize the time vector
             %Convert the time limits of the trace in datetime format
             Time1=datetime(startTime,'ConvertFrom','datenum');
             Time2=datetime(endTime,'ConvertFrom','datenum');
             %Compute the time vector of the  zlc analysis
             w=w_zlc;
             Time_zlc_1=Time1+seconds(w);
             Time_zlc_2=Time2;
             Time_zlc=Time_zlc_1:seconds(w):Time_zlc_2;
             time_zlc=[time_zlc;Time_zlc'];    
end

%==========================================================================
%%%%Plot of the Cross correlation coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting1(ii,ax)
    try
        nbins=str2num(get(h.nbins,'string'));
        timestep=hours(nbins);%time bin in hours
        w_zlc=str2num(get(h.w_zlc,'string'));%the analysis window 
        thr=str2num(get(h.thr,'string'));%the threshold value of the cross-correlation 
        cla(ax)%Clear the cross correlation axes
        Xcor=cell2mat(XCORR);
        Time=cell2mat(TIME_ZLC);
        Xcor=Xcor(:,ii);%Select the ii-th block of the cross correlation
        Time=Time(:,ii);%Select the ii-th block of the time
        idx=cell2mat(XCORR)>thr;%Selection of the values over the threshold
        idx=idx(:,ii);%Select the ii-th block of the selection
        Xcor=Xcor(idx);%Cross correlation selection of the values over the threshold
        Time=Time(idx);%Temporal selection of the values over the threshold
        %Set the the limits and steps analysis
        time1=min(Time);time2=max(Time)+datenum(timestep);xcor1=0;xcor2=1;xcorstep=0.1;
        %Calculation and plot the 3D temporal histograms of the Cross correlation coefficients
        DistrHist3D(ax,Time,Xcor,time1,time2,timestep,xcor1,xcor2,xcorstep)
        ylim(ax,[0 1])
        xlim(ax,[datetime(time1,'ConvertFrom','datenum') datetime(max(Time),'ConvertFrom','datenum')])%Convert x-axis limits in datetime format
        ylabel(ax,'Correlation');
        legend(ax,[num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
    catch
        %Error assessment
        if size(Xcor,1)~=size(Time,1);set(text1, 'String','Error. Incoherent dimensions of the ZLC array.');drawnow;return;end
    end
end

%==========================================================================
%%%%Plot of the BackAzimuth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting2(ii,ax)
    try
        nbins=str2num(get(h.nbins,'string'));
        w_zlc=str2num(get(h.w_zlc,'string'));%the analysis window  
        timestep=hours(nbins);%time bin in hours
        thr=str2num(get(h.thr,'string'));%the threshold value of cross-correlation 
        cla(ax)%Clear the backazimuth axes
        Az=cell2mat(Azimuth);
        Time=cell2mat(TIME_ZLC);
        Az=Az(:,ii);%Select the ii-th block of the backazimuth
        Time=Time(:,ii);%Select the ii-th block of the time
        idx=cell2mat(XCORR)>thr;%Selection of the values over the threshold
        idx=idx(:,ii);%Select the ii-th block of the selection
        Az=Az(idx);%Backazimuth selection of the values over the threshold
        Time=Time(idx);%Temporal selection of the values over the threshold
        %Set the the limits and steps analysis
        time1=min(Time);time2=max(Time)+datenum(timestep);az1=0;az2=360;azstep=10;
        %Calculation and plot the 3D temporal histograms of the BackAzimuth
        DistrHist3D(ax,Time,Az,time1,time2,timestep,az1,az2,azstep)
        ylim(ax,[0 360])
        xlim(ax,[datetime(time1,'ConvertFrom','datenum') datetime(max(Time),'ConvertFrom','datenum')])%Convert x-axis limits in datetime format
        ylabel(ax,'Back Azimuth (degrees)');
        legend(ax,[num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
    catch
        %Error assessment
        if size(Az,1)~=size(Time,1);set(text1, 'String','Error. Incoherent dimensions of the ZLC array.');drawnow;return;end
    end
end

%==========================================================================
%%%%Plot of the Ray parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting3(ii,ax)
    try
        nbins=str2num(get(h.nbins,'string'));
        w_zlc=str2num(get(h.w_zlc,'string'));%the analysis window   
        timestep=hours(nbins);%time bin in hours
        thr=str2num(get(h.thr,'string'));%the threshold value of cross-correlation 
        cla(ax)%Clear the ray parameter axes
        Ray=cell2mat(RayP);
        Time=cell2mat(TIME_ZLC);
        Ray=Ray(:,ii);%Select the ii-th block of the ray parameter 
        Time=Time(:,ii);%Select the ii-th block of the time
        idx=cell2mat(XCORR)>thr;%Selection of the values over the threshold
        idx=idx(:,ii);%Select the ii-th block of the selection
        Ray=Ray(idx);%Ray parameter selection of the values over the threshold
        Time=Time(idx);%Temporal selection of the values over the threshold
        %Set the the limits and steps analysis
        time1=min(Time);time2=max(Time)+datenum(timestep);ray1=0;ray2=10;raystep=0.1;
        %Calculation and plot the 3D temporal histograms of the ray parameter
        DistrHist3D(ax,Time,Ray,time1,time2,timestep,ray1,ray2,raystep)
        ylim(ax,[0 5])
        xlim(ax,[datetime(time1,'ConvertFrom','datenum') datetime(max(Time),'ConvertFrom','datenum')])%Convert x-axis limits in datetime format
        ylabel(ax,'RayParameter(s/km)');
        legend(ax,[num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
    catch
        %Error assessment
        if size(Ray,1)~=size(Time,1);set(text1, 'String','Error. Incoherent dimensions of the ZLC array.');drawnow;return;end
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
%%%%Set the Error analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function errormode(h, ~)
     text1.String=[];set(text1);drawnow;%Clear the Window message
    if (get(h,'Value') == get(h,'Max'))
        %Command message
        set(text1, 'String', 'Error mode is activated.');
        drawnow
        %Active Error analysis
        emod='y';
    else
        %Command message
        set(text1, 'String', 'Error mode is deactivated.');
        drawnow
        %Deactive Error analysis
        emod='n';
    end
end

%==========================================================================
%%%%Refresh the ZLC axes through the pre button%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loca1(~,~)
   %Error control
   if isempty(cell2mat(Azimuth));set(text1, 'String','Error. ZLC matrix is empty.');return;end
   %Decrease the slidval value
   slidval=slidval-1;if slidval<=0;slidval=1;end
   %Plot the new block of cross correlation
   ploting1(slidval,axes1);
   %Plot the new block of backazimuth
   ploting2(slidval,axes2);
   %Plot the new block of ray parameter
   ploting3(slidval,axes3);
        
end

%==========================================================================
%%%%Refresh the ZLC axes through the next button%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loca2(~,~)
    %Error control
    if isempty(cell2mat(Azimuth));set(text1, 'String', 'Error. ZLC matrix is empty.');return;end
    %Increase the slidval value
    slidval=slidval+1;if slidval>length(F)-1;slidval=length(F)-1;end
    %Plot the new block of cross correlation
    ploting1(slidval,axes1);
    %Plot the new block of backazimuth
    ploting2(slidval,axes2);
    %Plot the new block of ray parameter
    ploting3(slidval,axes3);
        
end

%==========================================================================
%%%%Set the limits of the y-axis on the cross correlation axes%%%%%%%%%%%%%
function setlimits1(~,~)
    try
        %Manual settings
        %Acquisition limits
        cor1_lim=str2num(get(h.cor1_lim,'string'));
        cor2_lim=str2num(get(h.cor2_lim,'string')); 
        ylim(axes1,[cor1_lim cor2_lim])
    catch
        %Default settings
        ylim(axes1,[0 1])
        %Reset cross correlation limits
        set(h.cor1_lim,'string','0');
        set(h.cor2_lim,'string','1');
    end
end

%==========================================================================
%%%%Set the limits of the y-axis on the backazimuth axes%%%%%%%%%%%%%%%%%%%
function setlimits2(~,~)
    try
        %Manual settings
        %Acquisition degrees limits
        az1_lim=str2num(get(h.az1_lim,'string'));
        az2_lim=str2num(get(h.az2_lim,'string')); 
        ylim(axes2,[az1_lim az2_lim])
    catch 
        %Default settings
        ylim(axes2,[0 360])
        %Reset backazimuth limits
        set(h.az1_lim,'string','0');
        set(h.az2_lim,'string','360');
    end
end

%==========================================================================
%%%%Set the limits of the y-axis on the ray parameter axes%%%%%%%%%%%%%%%%%
function setlimits3(~,~)
    try
        %Manual settings
        %Acquisition slowness limits
        ray1_lim=str2num(get(h.ray1_lim,'string'));
        ray2_lim=str2num(get(h.ray2_lim,'string')); 
        ylim(axes3,[ray1_lim ray2_lim])
    catch
        %Default settings
        ylim(axes3,[0 5])
        %Reset ray parameters limits
        set(h.ray1_lim,'string','0');
        set(h.ray2_lim,'string','5');
    end
end

%==========================================================================
%%%%Set the limits of the x-axis on the ZLC axes%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setlimits4(~,~)
    try
        %Manual settings
        %Acquisition time limits
        t1_lim=datetime(get(h.t1_lim,'string'),'InputFormat','yyMMdd-HHmmSS');%Convert in datetime format
        t2_lim=datetime(get(h.t2_lim,'string'),'InputFormat','yyMMdd-HHmmSS'); %Convert in datetime format
        xlim(axes1,[t1_lim t2_lim])
        xlim(axes2,[t1_lim t2_lim])
        xlim(axes3,[t1_lim t2_lim])
    catch
        %Default settings
        if ~isempty(cell2mat(TIME_ZLC)) 
            ii=slidval;
            Time=datetime(cell2mat(TIME_ZLC),'ConvertFrom','datenum');%Convert in datetime format
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
    w_zlc=str2num(get(h.w_zlc,'string'));   
    f1=str2num(get(h.f1,'string'));
    f2=str2num(get(h.f2,'string'));
    step_zlc=str2num(get(h.step_zlc,'string'));
    vel_zlc=str2num(get(h.vel_zlc,'string'));
    nbins=str2num(get(h.nbins,'string'));
    thr=str2num(get(h.thr,'string'));

    %Control and eventually reset parameters
    if (size(w_zlc,1)==0);set(h.w_zlc,'string',params.w_zlc);end
    if (size(f1,1)==0);set(h.f1,'string',params.f1_zlc);end
    if (size(f2,1)==0);set(h.f2,'string',params.f2_zlc);end
    if (size(step_zlc,1)==0);set(h.step_zlc,'string',params.step_zlc);end
    if (size(vel_zlc,1)==0);set(h.vel_zlc,'string',params.vel_zlc);end
    if (size(nbins,1)==0);set(h.nbins,'string','3');end
    if (size(thr,1)==0);set(h.thr,'string','0.75');end
end

%==========================================================================
%%%%Save ZLC plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save1(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    try
        nbins=str2num(get(h.nbins,'string'));
        w_zlc=str2num(get(h.w_zlc,'string'));%the analysis window  
        timestep=hours(nbins);%time bin in hours
        thr=str2num(get(h.thr,'string'));%the threshold value of the cross-correlation 
        ii=slidval;
        Az=cell2mat(Azimuth);
        Xcor=cell2mat(XCORR);
        Ray=cell2mat(RayP);
        Time=cell2mat(TIME_ZLC);
        Az=Az(:,ii);%Select the ii-th block of the backazimuth
        Xcor=Xcor(:,ii);%Select the ii-th block of the cross correlation
        Ray=Ray(:,ii);%Select the ii-th block of the ray parameter
        Time=Time(:,ii);%Select the ii-th block of the time
        idx=cell2mat(XCORR)>thr;%Selection of the values over the threshold
        idx=idx(:,ii);%Select the ii-th block of the selection
        Az=Az(idx);%Backazimuth selection of the values over the threshold
        Xcor=Xcor(idx);%Cross correlation selection of the values over the threshold
        Ray=Ray(idx);%Ray parameter selection of the values over the threshold
        Time=Time(idx);%Temporal selection of the values over the threshold
        %Set the the limits and steps analysis
        time1=min(Time);time2=max(Time)+datenum(timestep);xcor1=0;xcor2=1;xcorstep=0.1;
        az1=0;az2=360;azstep=10;ray1=0;ray2=10;raystep=0.1;
        new_figure=figure('Name','ZLC plot','NumberTitle','off');
        %Plot the values of the cross correlation
        subplot(3,1,1)
        %Calculation and plot the 3D temporal histograms of the Cross correlation coefficients
        NDistrHist3D(Time,Xcor,time1,time2,timestep,xcor1,xcor2,xcorstep)
        ylim([0 1])
        xlim([datetime(time1,'ConvertFrom','datenum') datetime(max(Time),'ConvertFrom','datenum')])%Convert x-axis limits in datetime format
        ylabel('Correlation');
        legend([num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
        %Plot the values of the backazimuth
        subplot(3,1,2)
        %Calculation and plot the 3D temporal histograms of the Backazimuth
        NDistrHist3D(Time,Az,time1,time2,timestep,az1,az2,azstep)
        ylim([0 360])
        xlim([datetime(time1,'ConvertFrom','datenum') datetime(max(Time),'ConvertFrom','datenum')])%Convert x-axis limits in datetime format
        ylabel('Back Azimuth (degress)');
        legend([num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
        %Plot the values of the ray parameter
        subplot(3,1,3)
        %Calculation and plot the 3D temporal histograms of the Ray parameter
        NDistrHist3D(Time,Ray,time1,time2,timestep,ray1,ray2,raystep)
        ylim([0 5])
        xlim([datetime(time1,'ConvertFrom','datenum') datetime(max(Time),'ConvertFrom','datenum')])%Convert x-axis limits in datetime format
        ylabel('RayParameter (s/km)');
        legend([num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])  
    catch
        %Error assessment
        if isempty(Time) | isempty (Az);set(text1, 'String','Error. ZLC matrix is empty.');drawnow;return;end
        if size(Az,1)~=size(Time,1);set(text1, 'String','Error. Incoherent dimensions of the ZLC array.');drawnow;return;end
    end
end

%==========================================================================
%%%%Plot the 3D Multivariate analysis Ray parameter vs. Backazimuth%%%%%%%%
function save2(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    try
        thr=str2num(get(h.thr,'string'));%the threshold value of the cross-correlation 
        ii=slidval;
        Az=cell2mat(Azimuth);
        Ray=cell2mat(RayP);
        Az=Az(:,ii);%Select the ii-th block of the backazimuth
        Ray=Ray(:,ii);%Select the ii-th block of the ray parameter
        idx=cell2mat(XCORR)>thr;%Selection of the values over the threshold
        idx=idx(:,ii);%Select the ii-th block of the selection
        Az=Az(idx);%Backazimuth selection of the values over the threshold
        Ray=Ray(idx);%Ray parameter selection of the values over the threshold
        %Set the the limits and steps analysis
        az1=0;az2=360;azstep=10;ray1=0;ray2=10;raystep=0.1;
        new_figure=figure('Name','Ray Parameter vs. BackAzimuth','NumberTitle','off');
        %Calculation and plot the 3D Multivariate analysis
        NewHist3D(Ray,Az,ray1,ray2,raystep,az1,az2,azstep)
        set(gca,'ytick', [0:36:360], 'ylim', [0 360])
        set(gca,'xtick', [0:0.5:5], 'xlim', [0 5])
        ylabel('Backazimuth (degrees)');
        xlabel('RayParameter (s/km)');
        legend([num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz']) 
    catch
        %Error assessment
        if isempty(Ray) | isempty (Az);set(text1, 'String','Error. ZLC matrix is empty.');drawnow;return;end
        if size(Az,1)~=size(Ray,1);set(text1, 'String','Error. Incoherent dimensions of the ZLC array.');drawnow;return;end
    end
end

%==========================================================================
%%%%Plot the 3D Multivariate analysis of the horizontal slowness%%%%%%%%%%%
function save3(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    try
        thr=str2num(get(h.thr,'string'));%the threshold value of the cross-correlation 
        ii=slidval;
        sx=cell2mat(ssx);
        sy=cell2mat(ssy);
        sx=sx(:,ii);sy=sy(:,ii);%Select the ii-th block of the horizontal slowness
        idx=cell2mat(XCORR)>thr;%Selection of the values over the threshold
        idx=idx(:,ii);%Select the ii-th block of the selection
        sx=sx(idx);sy=sy(idx);%Slowness selection of the values over the threshold
        new_figure=figure('Name','Slowness Grid','NumberTitle','off');
        %Calculation and plot the 3D Multivariate analysis
        NewHist3D(sx,sy,-5,5,0.1,-5,5,0.1)
        set(gca,'ytick', [-5:1:5], 'ylim', [-5 5])
        set(gca,'xtick', [-5:1:5], 'xlim', [-5 5])
        ylabel('Sy (s/km)');
        xlabel('Sx (s/km)');
        legend([num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz']) 
    catch
        %Error assessment
        if isempty(sx) | isempty (sy);set(text1, 'String','Error. ZLC matrix is empty.');drawnow;return;end
        if size(sx,1)~=size(sy,1);set(text1, 'String','Error. Incoherent dimensions of the ZLC array.');drawnow;return;end
    end
end

%==========================================================================
%%%%Plot the statistics of the Incidence and Backazimuth%%%%%%%%%%%%%%%%%%%
function save4(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    try
        %Activation/Deactivation of the error boxplots
        if emod=='y' & ~isempty(cell2mat(Eb))
            thr=str2num(get(h.thr,'string'));%the threshold value of the cross-correlation 
            ii=slidval;
            Inc=cell2mat(Incidence);
            Az=cell2mat(Azimuth);
            Einc=cell2mat(Ei);
            Ebaz=cell2mat(Eb);
            Inc=Inc(:,ii);%Select the ii-th block of the incidence
            Az=Az(:,ii);%Select the ii-th block of the backazimuth
            Einc=Einc(:,ii);%Select the ii-th block of the incidence error
            Ebaz=Ebaz(:,ii);%Select the ii-th block of the backazimuth error
            idx=cell2mat(XCORR)>thr;%Selection of the values over the threshold
            idx=idx(:,ii);%Select the ii-th block of the selection
            Inc=Inc(idx);%Incidence selection of the values over the threshold
            Az=Az(idx);%Azimuth selection of the values over the threshold
            Einc=Einc(idx);%Error incidence selection of the values over the threshold
            Ebaz=Ebaz(idx);%Error backazimuth selection of the values over the threshold
            %Incidence bins
            bins=0:10:90;
            %Backazimuth bins
            ros=Az;
            fin=(2*pi*ros)./360;
            bins2=[0:5*((2*pi)/360):355*((2*pi)/360)];
            new_figure=figure('Name','Incidence vs. BackAzimuth','NumberTitle','off');
            %Plot the histogram of the incidence
            subplot(2,2,1)
            newh=histogram(Inc,bins,'FaceColor','b');
            newh.Normalization = 'probability';
            xlim([0 90])
            ylabel('Probability')
            legend([num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
            grid on; grid minor;
            title(['Incidence (degrees)'])
            %Plot the polar histogram of the backazimuth
            subplot(2,2,2)
            newh2=polarhistogram(fin,bins2);
            newh2.Normalization = 'probability';
            set(gca,'ThetaZeroLocation','top',...
            'ThetaDir','clockwise');
            title(['Backazimuth (degrees)'])
            legend([num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
            %Plot the distribution of the incidence error
            subplot(2,2,3)
            boxplot(Einc)
            grid on; grid minor;
            ylabel('Error Incidence (degrees)')
            set(gca,'xtick',[])
            %Plot the distribution of the backazimuth error
            subplot(2,2,4)
            boxplot(Ebaz)
            grid on; grid minor;
            ylabel('Error Back Azimuth (degrees)')
            set(gca,'xtick',[])
        else
            thr=str2num(get(h.thr,'string'));%the threshold value of the cross-correlation 
            ii=slidval;
            Inc=cell2mat(Incidence);
            Az=cell2mat(Azimuth);
            Inc=Inc(:,ii);%Select the ii-th block of the incidence
            Az=Az(:,ii);%Select the ii-th block of the backazimuth
            idx=cell2mat(XCORR)>thr;%Selection of the values over the threshold
            idx=idx(:,ii);%Select the ii-th block of the selection
            Inc=Inc(idx);%Incidence selection of the values over the threshold
            Az=Az(idx);%Backzimuth selection of the values over the threshold
            %Incidence bins
            bins=0:10:90;
            %Backazimuth bins
            ros=Az;
            fin=(2*pi*ros)./360;
            bins2=[0:5*((2*pi)/360):355*((2*pi)/360)];
            new_figure=figure('Name','Incidence vs. BackAzimuth','NumberTitle','off');
            %Plot the histogram of the incidence
            subplot(2,1,1)
            newh=histogram(Inc,bins,'FaceColor','b');
            newh.Normalization = 'probability';
            xlim([0 90])
            ylabel('Probability')
            legend([num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
            grid on; grid minor;
            title(['Incidence (degrees)'])
            %Plot the polar histogram of the backazimuth
            subplot(2,1,2)
            newh2=polarhistogram(fin,bins2);
            newh2.Normalization = 'probability';
            set(gca,'ThetaZeroLocation','top',...
            'ThetaDir','clockwise');
            title(['Backazimuth (degrees)'])
            legend([num2str(F(ii)) '-' num2str(F(ii+1)) ' Hz'])
        end
    catch
        %Error assessment
        if isempty(Az);set(text1, 'String','Error. ZLC matrix is empty.');drawnow;return;end
        if isempty(cell2mat(Eb)) | isempty (cell2mat(Ei));set(text1, 'String','Error. Error matrix is empty.');drawnow;return;end
        if size(Az,1)~=size(Ray,1);set(text1, 'String','Error. Incoherent dimensions of the ZLC array.');drawnow;return;end   
    end
end

%==========================================================================
%%%%Save the results of the zero lag cross correlation analysis%%%%%%%%%%%%
function save_zlc(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Validation of the results
    if isempty(params.save);set(text1, 'String', sprintf('Error. No output directory.'));drawnow;return;end  
    if isempty(cell2mat(Azimuth));set(text1, 'String', sprintf('Error. ZLC matrix is empty.'));drawnow;return;end
    if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
    %Loop through the frequency analysis
    for jj=1:size(Azimuth,2)
        %Loop through the number of Input files (ff)
        for kk=1:size(Azimuth,1)
            time=TIME_ZLC{kk,jj};baz=Azimuth{kk,jj};rayp=RayP{kk,jj};xcmax=XCORR{kk,jj};sx=ssx{kk,jj};sy=ssy{kk,jj};Ebaz=Eb{kk,jj};Eray=Er{kk,jj};Einc=Ei{kk,jj};inc=Incidence{kk,jj};
            if isempty(Ebaz);Ebaz=NaN(length(baz),1);Eray=NaN(length(baz),1);Einc=NaN(length(baz),1);end
            freq1=F(jj);freq2=F(jj+1);Ffolder=[num2str(freq1) '-' num2str(freq2) 'Hz'];
            name=ff{kk};[~,name,~] = fileparts(name);
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
            pathout=[params.save '/ZLC/' Ffolder '/' jday '/'];if ~exist(pathout, 'dir');mkdir(pathout);end 
            output=[pathout name params.outputformat];
            %Create a table of the results
            data=table(time,baz,rayp,inc,Ebaz,Eray,Einc,sx,sy,xcmax);data.Properties.VariableNames={'Time','BackAzimuth', 'RayParameter','Incidence','ErrorBaz','ErrorRay','ErrorInc','Sx','Sy','Xcoef'};
            data.Properties.VariableUnits = {'datenum','degrees','s/km ','degrees','degrees','s/km ','degrees','s/km ','s/km ',' '};
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
    set(text1, 'String', sprintf(['ZLC' params.outputformat '-Observed data are saved.']));
    drawnow
end
end