%%%%%%%%%%%%% Calculation of the Radial Semblance function %%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params: user parameters
% avfig: Radial Semblance figure
% panel_color: colors of the main panel 
% fig_color: colors of the figure
% pan1: control right side panel 
% text1: Window message
% axes1,axes2,axes3: axes of the current figure
% x1_lim, x2_lim: limit of x-axis on XYZ axes
% y1_lim, y2_lim: limit of y-axis on XYZ axes
% z1_lim, z2_lim: limit of z-axis on XYZ axes
% slidval: value of dynamic plots 
% h: UI control parameters
%%%%%%%%%%%h.variables/params.variables/variables%%%%%%%%%%%%%%%%%%%%%%%%%%
% rad: radio button within a button group 
% rad1, rad2, rad3, rad4: button groups 
%-------------------------------- General settings
%--- General info
% sta_ref: the name of the reference station 
% component, comp: the name of the component/channel
% comp_chan: the station system used (component/channel)
% station: station names
% NStation: number of stations used
% coord: the path of the coordinate file (.mat) 
% coord_sys: the coordinate system (metrics/degrees)
% xs/XS, ys/YS, zs/ZS, utm_zone: coordinates in UTM system
% LAT, LON, ELE: coordinates in wgs84 system
%--- Data file
% data: the Input folder
% save: the Output folder
% outputformat: the format (.mat/.txt) of the output files
% search_mat: term of research file 
% pathname/pathnam/pathaname2, fname: Input folder and list of filenames
% ff: list of files with their path
% list: list of three components directories
% list2: list of seismic traces directories
% mytrace: seismic trace structure
% dt: sampling interval
% fs: sampling rate
% nyq: nyquist number
% z: amplitude vector of the seismic trace
% nt: number of samples
% startTime, endTime: temporal limits of the traces 
%--- Semblance setting
% delim: the type of analysis (Automatic continous/Automatic discrete/Manual discrete)
% emod: Error analysis routine (on/off)
% w_sem: the analysis window for the Radial Semblance analysis (s)
% wl: analysis window (samples)
% vel_sem: the velocity value (km/s)
% step_sem: the grid step size (km)
% topo_corr: the topographic correction (yes o no)
% freq_sem: the  frequency range for the band-pass filter (Hz)
% xmin: the minimum X-limit of the grid (km)
% xmax: the maximum X-limit of the grid (km)
% ymin: the minimum Y-limit of the grid (km)
% ymax: the maximum Y-limit of the grid (km)
% zmin: the minimum Z-limit of the grid (km)
% zmax: the maximum Z-limit of the grid (km)
% valor: maximum Semblance value
% ix1,ix2,ix3: indexes of the maximum Radial Semblance value
% staidx: station index of reference
% XXq/Xq,YYq/Yq,ZZq/Zq: 3D search grid Nx*Ny*Nz
% xxq/xq,yyq/yq,zzq/zq: coordinate vectors of the search grid
% LOC,LOCAL: index matrix of the Localization vectors
% Err,ERROR: Analysis errors matrixs (km)
% loc,GRID: Semblance grid matrix
% x,y,z: Localization vectors(km)
% Ex,Ey,Ez: Average values of the localization vectors(km) 
% Vex, VEy, VEz: Average values of the errors matrix(km)
% vv: dimension of the scatters
% Name: filename vector
% ttrig,trig,TRIGG: picking times
% ttit: initial title
% tit: final title
% AA/A,XX/X,YYy/Y: the DEM variables
% J_error: JackKnife results (km)
% Z/ZZ,N/NN,E/EE: Three components matrixs

function RadSemLoc()
%% Declaration global parameters
clc
global params Xq Yq Zq xq yq zq AA XX YYy fs z XXq YYq ZZq xxq yyq zzq ttit
comp=params.component;
component=[];
ff= cell(1);
slidval=1;
pathname=[];pathname2=[];pathnam=[];
LOCAL=[];
ERROR=[];
GRID=[];
Name=[];
TRIGG=[];
tit=[];
LAT=[];LON=[];ELE=[];xs=[];ys=[];zs=[];station={}; utm_zone=[];emod='n';
%% Radial Semblance figure
fig_color= [.94 .94 .94];
panel_color= [.97 .97 .97];

avfig= figure('Name', ' Radial Semblance Analysis', 'position',get(0,'screensize'), ...
    'NumberTitle','off','units','normalized', 'color', fig_color,'MenuBar','none');

%Window message to user
text1= uicontrol('Style','edit','FontName','Arial','Fontsize',10,'Fontweight','bold',...
                'String','Window Message','Max',5,'Min',0,'HorizontalAlignment','left',...
                'BackgroundColor', 'w','Enable','inactive',...
                'units','normalized','Position',[.05 .03 .3 .10]);


%% Axeses
% XY plane axes
axes1= axes('position',[.05 .25 .32 .6]);
hold(axes1,'on')
grid(axes1,'on')
grid(axes1,'minor')
colorbar(axes1)
colormap(axes1,'parula')
%Set axis labels and title text on the basis of the type of coordinates system (UTM/wgs84 system)
switch params.coord_sys
    case 'm'
        xlabel(axes1,'X (m)')
        ylabel(axes1,'Y (m)')
        ttit1=cellstr('yyyy.jjj.hhmmss');
        ttit2=cellstr(sprintf('Max point (%.3f) : Nx=%.2f, Ny=%.2f, Nz=%.3f \n',0,0,0,0));
        ttit3=cellstr(sprintf('Error Max point  : Ex=%.3f, Ey=%.3f, Ez=%.3f \n',0,0,0));
        ttit=[ttit1;ttit2;ttit3];
    case 'degrees'
        ylabel(axes1,'Latitude (decimal degrees)')
        xlabel(axes1,'Longitude (decimal degrees)')
        ttit1=cellstr('yyyy.jjj.hhmmss');
        ttit2=cellstr(sprintf('Max point (%.3f) : Lat=%.3f, Lon=%.3f, Elv=%.3f \n',0,0,0,0));
        ttit3=cellstr(sprintf('Error Max point  : Ex=%.4f, Ey=%.4f, Ez=%.3f \n',0,0,0));
        ttit=[ttit1;ttit2;ttit3];
end
title(axes1,ttit);

% YZ plane axes
axes2= axes('position',[0.4 0.68 0.4 0.28]);
hold(axes2,'on')
grid(axes2,'on')
grid(axes2,'minor')
colorbar(axes2)
colormap(axes2,'parula')
%Set axis labels on the basis of the type of coordinates system (UTM/wgs84 system)
switch params.coord_sys
    case 'm'
        xlabel(axes2,'Y (m)')
        ylabel(axes2,'Z (m)')
    case 'degrees'
        xlabel(axes2,'Latitude (decimal degrees)')
        ylabel(axes2,'Elevation (km)') 
end

% XZ plane axes
axes3= axes('position',[0.4 0.14 0.4 0.28]);
hold(axes3,'on')
grid(axes3,'on')
grid(axes3,'minor')
colorbar(axes3)
colormap(axes3,'parula')
%Set axis labels on the basis of the type of coordinates system (UTM/wgs84 system)
switch params.coord_sys
    case 'm'
        xlabel(axes3,'X (m)')
        ylabel(axes3,'Z (m)')
    case 'degrees'
        xlabel(axes3,'Longitude(decimal degrees)')
        ylabel(axes3,'Elevation (km)')       
end

%% Control right side panel
pan1= uipanel(avfig,'visible','on','Position',[.81 .03 .185 .94],...
    'BackgroundColor', panel_color);
uipanel(avfig,'visible','on','Position',[.81 .27 .185 .01],...
    'BackgroundColor', fig_color);

%% Text boxes
%-------------------------------- Radial Semblance settings
uicontrol(pan1,'Style','text', 'String','RadSemb setting',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.03 .76 .6 .06],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Win (s)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .75 .4 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Vel (km/s)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .68 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Freq (Hz)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .63 .4 .04],...
    'BackgroundColor',panel_color);

%-------------------------------- Grid settings
uicontrol(pan1,'Style','text', 'String','GRID setting',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.55 .76 .6 .06],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Xmin (km)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .75 .4 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Xmax (km)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .70 .4 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Ymin (km)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .65 .4 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Ymax (km)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .58 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Zmin (km)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .53 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Zmax (km)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .48 .4 .04],...
    'BackgroundColor',panel_color);

%-------------------------------- Axis settings
uicontrol(pan1,'Style','text', 'String','Axis setting',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.03 .42 .4 .04],...
    'BackgroundColor',panel_color);
%%Create the text boxes on the basis of the type of coordinates system (UTM/wgs84 system)
switch params.coord_sys
    case 'degrees'
        uicontrol(pan1,'Style','text', 'String','Lon1 (deg)',...
            'HorizontalAlignment','left',...
            'Units','normalized','Position',[.02 .39 .4 .02],...
            'BackgroundColor',panel_color);
        uicontrol(pan1,'Style','text', 'String','Lat1 (deg)',...
            'HorizontalAlignment','left',...
            'Units','normalized','Position',[.02 .34 .4 .02],...
            'BackgroundColor',panel_color);
        uicontrol(pan1,'Style','text', 'String','Elv1 (deg)',...
            'HorizontalAlignment','left',...
            'Units','normalized','Position',[.02 .27 .4 .04],...
            'BackgroundColor',panel_color);
        uicontrol(pan1,'Style','text', 'String','Lon2 (deg)',...
            'HorizontalAlignment','left',...
            'Units','normalized','Position',[.54 .37 .4 .04],...
            'BackgroundColor',panel_color);
        uicontrol(pan1,'Style','text', 'String','Lat2 (deg)',...
            'HorizontalAlignment','left',...
            'Units','normalized','Position',[.54 .32 .4 .04],...
            'BackgroundColor',panel_color);
        uicontrol(pan1,'Style','text', 'String','Elv2 (deg)',...
            'HorizontalAlignment','left',...
            'Units','normalized','Position',[.54 .27 .4 .04],...
            'BackgroundColor',panel_color);
    case 'm'
        uicontrol(pan1,'Style','text', 'String','x1 (km)',...
            'HorizontalAlignment','left',...
            'Units','normalized','Position',[.03 .39 .4 .02],...
            'BackgroundColor',panel_color);
        uicontrol(pan1,'Style','text', 'String','y1 (km)',...
            'HorizontalAlignment','left',...
            'Units','normalized','Position',[.03 .34 .4 .02],...
            'BackgroundColor',panel_color);
        uicontrol(pan1,'Style','text', 'String','z1 (km)',...
            'HorizontalAlignment','left',...
            'Units','normalized','Position',[.03 .27 .4 .04],...
            'BackgroundColor',panel_color);
        uicontrol(pan1,'Style','text', 'String','x2 (km)',...
            'HorizontalAlignment','left',...
            'Units','normalized','Position',[.55 .37 .4 .04],...
            'BackgroundColor',panel_color);
        uicontrol(pan1,'Style','text', 'String','y2 (km)',...
            'HorizontalAlignment','left',...
            'Units','normalized','Position',[.55 .32 .4 .04],...
            'BackgroundColor',panel_color);
        uicontrol(pan1,'Style','text', 'String','z2 (km)',...
            'HorizontalAlignment','left',...
            'Units','normalized','Position',[.55 .27 .4 .04],...
            'BackgroundColor',panel_color);   
end

%% Editable text 
%-------------------------------- Radial Semblance settings
h.w_sem= uicontrol(pan1,'Style','edit', 'String',params.w_sem,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .74 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.vel_sem= uicontrol(pan1,'Style','edit', 'String',params.vel_sem,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .69 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.freq_sem= uicontrol(pan1,'Style','edit', 'String',params.freq_sem,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .64 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.delim= uicontrol(pan1,'Style','popup', ...
    'String',{'Continuous','Auto-Discrete','Man-Discrete'},'Value',2,...
    'HorizontalAlignment','left','Enable','on',...
    'Units','normalized','Position',[.58 .81 .4 .1],...
    'BackgroundColor','w');

%-------------------------------- Grid settings
h.xmin= uicontrol(pan1,'Style','edit', 'String',params.xmin,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .74 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.xmax= uicontrol(pan1,'Style','edit', 'String',params.xmax,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .69 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.ymin= uicontrol(pan1,'Style','edit', 'String',params.ymin,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .64 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.ymax= uicontrol(pan1,'Style','edit', 'String',params.ymax,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .59 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.zmin= uicontrol(pan1,'Style','edit', 'String',params.zmin,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .54 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.zmax= uicontrol(pan1,'Style','edit', 'String',params.zmax,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .49 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);

%-------------------------------- Axis settings
h.x1_lim= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .38 .2 .04],...
    'BackgroundColor','w','callback',@setlimits1);
h.y1_lim= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .33 .2 .04],...
    'BackgroundColor','w','callback',@setlimits2);
h.z1_lim= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .28 .2 .04],...
    'BackgroundColor','w','callback',@setlimits3);
h.x2_lim= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .38 .2 .04],...
    'BackgroundColor','w','callback',@setlimits1);
h.y2_lim= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .33 .2 .04],...
    'BackgroundColor','w','callback',@setlimits2);
h.z2_lim= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .28 .2 .04],...
    'BackgroundColor','w','callback',@setlimits3);


%% Add Radio buttons
% Editable text
switch params.comp_chan
    case 'Comp'
        uicontrol(pan1,'Style','text', 'String','Component',...
            'HorizontalAlignment','left',...
            'Units','normalized','Position',[.05 .85 .4 .05],...
            'BackgroundColor',panel_color);
    case 'Chan'
        uicontrol(pan1,'Style','text', 'String','Channel',...
            'HorizontalAlignment','left',...
            'Units','normalized','Position',[.05 .85 .4 .05],...
            'BackgroundColor',panel_color);
end

% Component selection radio buttons
h.rad= uibuttongroup(pan1,'units','normalized','BackgroundColor',panel_color,...
    'bordertype','none','Position',[.08 .82 .5 .05]);

set(h.rad,'SelectionChangeFcn',@radcbk);

h.rad1 = uicontrol( h.rad, 'Style','Radio','String','Z',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.0 0 .25 1],'HandleVisibility','off');
h.rad2 = uicontrol( h.rad, 'Style','Radio','String','N',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.25 0 .25 1],'HandleVisibility','off');
h.rad3 = uicontrol( h.rad, 'Style','Radio','String','E',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.5 0 .25 1],'HandleVisibility','off');
h.rad4 = uicontrol( h.rad, 'Style','Radio','String','F',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.75 0 .25 1],'HandleVisibility','off');


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

% Error mode checkbox           
uicontrol('Style','checkbox','BackgroundColor', panel_color,...
                'String','Error Mode','units','normalized','Tag','errmod',...
                'Value',0,'Position',[.9175 .8 .06 .05],'callback',@errormode);

% Editable texts
uicontrol(pan1,'Style','text', 'String','RSemblance',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.03 .17 .4 .07],...
    'BackgroundColor',panel_color);

% Calculation and save button
uicontrol('style','pushbutton', 'string','Calculation','units','normalized',...
    'position',[.815 .18 .08 .0494], 'callback',@calc); 
uicontrol('style','pushbutton', 'string','Save Main Plot','units','normalized',...
    'position',[.815 .115 .08 .0494], 'callback',@save1); 
uicontrol('style','pushbutton', 'string','Save Data','units','normalized',...
    'position',[.815 .05 .08 .0494], 'callback',@save_sem);
uicontrol('style','pushbutton', 'string','Volume','units','normalized',...
    'position',[.905 .18 .08 .0494], 'callback',@save2);
uicontrol('style','pushbutton', 'string','Scatter locs','units','normalized',...
    'position',[.905 .115 .08 .0494], 'callback',@save3);

% Next and pre button
uicontrol('style','pushbutton', 'string','Next','units','normalized',...
    'position',[.77 .5 .03 .04], 'callback',@loca2);
uicontrol('style','pushbutton', 'string','Pre','units','normalized',...
    'position',[.77 .56 .03 .04], 'callback',@loca1);

%% Other functions (nested)
%%%%Return the list of files manually selected%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function openfnc(~,~)
     LOCAL=[];ERROR=[];GRID=[];Name=[];TRIGG=[];%Initialize the main vectors/matrixs
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);cla(axes2);cla(axes3)%Clear the axes
    title(axes1,ttit);%Clear the title
    tt= 0;pathnam=[];pathname2=[];
    ff=cell(1);%Initialize the list of files
    %Set the term of research on the basis of the station system (component/channel)
    switch params.comp_chan
        case 'Comp'
            search_mat=strcat(params.sta_ref,'.BH',comp,'*.mat');
        case 'Chan'
            %Error message
            textLabel = sprintf('Error. No three components files are found.');
            set(text1, 'String', textLabel);
            drawnow
            return
    end
    [fname, pathname] = uigetfile({search_mat;'*.*'},'File Selector',...
        strcat(cd,'/Data/output/'),'MultiSelect','on'); %Selection of the files    
     component=comp; 
    
    [~,nnn]= size(fname);
    if nnn==1 && fname==0
        %Error message
        set(text1, 'String', 'Error. No files selection made.');
        drawnow
    else
        fname= cellstr(fname);
        [~,np]= size(fname);%Number of files selected
        if np<1
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
     LOCAL=[];ERROR=[];GRID=[];Name=[];TRIGG=[];%Initialize the main vectors/matrixs
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);cla(axes2);cla(axes3)%Clear the axes
    title(axes1,ttit);%Clear the title
    pathname=[];
    pathnam= uigetdir(strcat(cd,'/Data/output/'));
    pathnam=[pathnam '/'];
    %Set the term of research on the basis of the station system (component/channel)
    switch params.comp_chan
        case 'Comp'
            search_mat=strcat(params.sta_ref,'.BH',comp,'*.mat');
        case 'Chan'
            %Error meassage
            textLabel = sprintf('Error. No three components files are found.');
            set(text1, 'String', textLabel);
            drawnow
            return  
    end
    
    if pathnam~=0
    % Get all files in the directory and subdirectories
    [ff,pathname2]= getAllFiles(pathnam,search_mat);
     ff=ff';pathname2=pathname2;
    
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
%%%%Calculation and plot of the Radial Semblance function%%%%%%%%%%%%%%%%%%
function calc(~,~) 
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);cla(axes2);cla(axes3)%Clear the axes
    title(axes1,ttit);%Clear the title
    YY={};%Initialize the events vector
    station=[];%Initialize the station names vector
    try
        %Command message
        textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
        set(text1, 'String', textLabel); 
        drawnow
        %Acquisition and validation of the analysis parameters
        w_sem=str2num(get(h.w_sem,'string'));if w_sem<=0;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
        v=str2num(get(h.vel_sem,'string'));if v<=0;set(text1, 'String','Error. Invalid velocity value.');drawnow;return;end
        freq_sem=deblank(get(h.freq_sem,'string'));
        freq=split(freq_sem,'-');
        f1=str2num(freq{1});f2=str2num(freq{2});delim=get(h.delim,'Value');
        if f1<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        if f2<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        xmin=str2num(get(h.xmin,'string'));
        xmax=str2num(get(h.xmax,'string'));
        ymin=str2num(get(h.ymin,'string'));
        ymax=str2num(get(h.ymax,'string'));
        zmin=str2num(get(h.zmin,'string'));
        zmax=str2num(get(h.zmax,'string'));
        if xmin<=0;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
        if xmax<=0 | xmax<=xmin;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
        if ymin<=0;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
        if ymax<=0 | ymax<=ymin;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
        if zmin<=0;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
        if zmax<=0 | zmax<=zmin;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
        if params.step_sem<=0;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
        if strcmp(params.topo_corr,'yes') & strcmp(params.topo_corr,'no');set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
        %Load station coordinate file
        xs=[];ys=[];zs=[];station={};LAT=[];LON=[];ELE=[];%Initialize coordinates vectors
        load(params.coord)
        %Error control
        if isempty(params.coord) | isempty(station);set(text1, 'String','Error. No coordinate file directory or No correct coordinate file.');drawnow;return;end
        %Execute one of two groups of routines (coordinates in wgs84/UTM)
        switch params.coord_sys
            case 'degrees'
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
                %Error control
                if isempty(xs);set(text1, 'String','Error. No correct coordinate file.');drawnow;return;end
        end
        %Reading .mat files
        [~, len]=size (ff);%Number of files selected
        LOCAL=[];%Initialize the Localization matrix
        ERROR=[];%Initialize the Error matrix
        GRID=[];%Initialize the Semblance matrix
        Name=[];%Initialize the filename matrix
        %Construction of  the Search grid and the DEM grid (UTM system)
        dem_file='mapfile_m.mat';
        [Xq,Yq,Zq,Nx,Ny,Nz,xq,yq,zq,AA,XX,YYy]=DEMconst(dem_file,params.topo_corr,xmin,xmax,ymin,ymax,zmin,zmax,params.step_sem);
        NStation=length(station);%Number of stations used
        %Define the station index of reference
        staidx=find(strcmp(station,params.sta_ref));
        %Picking routine 
        [TRIGG,ZZ,EE,NN,XSS,YSS,ZSS,Name,fs]=SelectSignal(len,NStation,f1,f2,staidx,delim);
        w1=round(w_sem*fs);
        %Loop through the number of files selected
        for i=1:len
            %Select the i-th components of the signals
            Z=ZZ{i};
            E=EE{i};
            N=NN{i};
            %Select the i-th station coordinates
            XS=XSS{i};YS=YSS{i};ZS=ZSS{i};
            ttrig=TRIGG{i};%Select the i-th picking times
            ntrig=length(ttrig);%number of picking times
            Err=[];%Initialize the Error matrix
            Grid=0;%Initialize the Semblance grid
            %Loop through the number of pickings
            for k=1:ntrig
                trig=ttrig(k);%Select the k-th picking time
                %Application of the Radial Semblance algorithm
                [valor,ix1,ix2,ix3,loc]=RadialS(Z,E,N,Xq,Yq,Zq,xq,yq,zq,NStation,v,Nx,Ny,Nz,XS,YS,ZS,fs,trig,staidx,w1);
                J_error=[NaN NaN NaN];%Initialize the standard error of the JackKnife estimator
                %Active/Deactive the the JackKnife routine 
                if emod=='y'
                    [J_error]=RadJackKnife(Z,E,N,Xq,Yq,Zq,xq,yq,zq,NStation,v,Nx,Ny,Nz,XS,YS,ZS,fs,trig,staidx,w1,ix1,ix2,ix3);
                end
                %Concatenate the results
                Err=[Err;J_error];
                Grid=Grid+loc; 
            end
            %Averaging the results
            J_error=mean(Err,1);loc=Grid./ntrig;
            %Set the new points of maximum probability
            [valor,ix]=max(loc(:),[],'omitnan');
            [ix1, ix2,ix3] = ind2sub(size(loc),ix);
            
            %Concatenate the results
            LOC=[ix1,ix2,ix3,valor,NStation];
            LOCAL=[LOCAL;LOC];
            ERROR=[ERROR;J_error];
            GRID=[GRID;{loc}]; 
        end
        %Set the block of results
        slidval=len;
        %Plot the results on the basis of the type of coordinates system (UTM/wgs84 system)
        switch params.coord_sys
            case 'm'
                %Plot the results in the XY plane axes (UTM system) 
                ploting1(slidval,axes1)
                %Plot the results in the YZ plane axes (UTM system)
                ploting2(slidval,axes2)
                %Plot the results in the XZ plane axes (UTM system) 
                ploting3(slidval,axes3)
            case 'degrees'
                %Construction of  the Search grid and the DEM grid (wgs84 system)
                dem_file='mapfile_deg.mat';
                [XXq,YYq,ZZq,xxq,yyq,zzq,~,~,~]=DEMconstDeg(dem_file,params.topo_corr,xmin,xmax,ymin,ymax,zmin,zmax,params.step_sem,utm_zone);
                %Plot the results in the XY plane axes (wgs84 system) 
                ploting1_deg(slidval,axes1)
                %Plot the results in the YZ plane axes (wgs84 system)
                ploting2_deg(slidval,axes2)
                %Plot the results in the XZ plane axes (wgs84 system) 
                ploting3_deg(slidval,axes3)
        end
        %Command window
        set(text1, 'String','Calculation finished.');
        drawnow
    catch
        %Error assessment 
        if isempty(params.coord) | isempty(station);set(text1, 'String','Error. No coordinate file directory or No correct coordinate file.');drawnow;return;end
        if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
        if isempty(z);set(text1, 'String','Error. No correct signals data.');drawnow;return;end
        if length(z)<=round(w_sem*fs);set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end
        if f1>=fs/2 | f2>fs/2 | f0>=fs/2;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        if isempty(LOCAL);set(text1, 'String','Error. Localization matrix is empty.');drawnow;return;end
        if size(LOCAL,1)~=len;set(text1, 'String','Error. Invalid picking routine.');drawnow;return;end
    end
end

%==========================================================================
%%%%%%%%%%%%%%%%Picking routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TRIG,Z,E,N,XSS,YSS,ZSS,Namemain,fs]=SelectSignal(len,NStation,f1,f2,staidx,delim)
    try
        Z=[];E=[];N=[];%Initialize the components matrixs
        TRIG=[];%Initialize the picking times vectors
        XSS=[];YSS=[];ZSS=[];%Initialize the station coordinates vectors
        mytrace=[];%Initialize the seismic trace structure
        Namemain=[];%Initialize the filename vectors
        w_sem=str2num(get(h.w_sem,'string'));%Get the analysis window
        %Loop through the number of files selected
        for i=1:len
            file_name=ff{i};%return the i-th file
            %Create the term of the signals research
            switch params.comp_chan 
                case 'Comp'
                    nname=file_name(end-18:end-4);
                case 'Chan'
                    %Error message
                     textLabel = sprintf('Error. No three components files are found.');
                     set(text1, 'String', textLabel);
                     drawnow
                     return  
            end
            %Create the term of the components research
            name=file_name(end-18:end);
            %Control and set the Input folder
            if ~isempty(pathnam);pathname=pathname2{i};pathname=[pathname '/'];end
            %Search all seismic traces based on the station names
            list2=searchAllStaz(pathname,station,component,nname,params.comp_chan);nlist2=length(list2);
            %Control the calculation through the number of stations
            if nlist2==NStation
                Zz=[];Ee=[];Nn=[];NameE=[];NameN=[];NameZ=[];%Initialize the main vectors
                %Loop through the number of stations used
                for jj=1:nlist2
                    %Search all three component per trace
                    list=dir(strcat(pathname,'/',list2(jj).name(1:4),'*',list2(jj).name(10:end)));
                    %Control the calculation through the number of components
                    if length(list)==3
                        %Load the three components 
                        fileE=strcat(pathname,'/',list(1).name);
                        fileN=strcat(pathname,'/',list(2).name);
                        fileZ=strcat(pathname,'/',list(3).name);
                        load(fileE)
                        EE=mytrace.data;
                        load(fileN)
                        NN=mytrace.data;
                        load(fileZ)
                        ZZ=mytrace.data;fs=mytrace.sampleRate;
                        dt=1/fs;fn=fs/2;
                        %Filter the seismic traces with a band-pass
                        [Bf,Af]=butter(2,[f1 f2]/(fs/2));
                        ZZ=filter(Bf,Af,ZZ);
                        EE=filter(Bf,Af,EE);
                        NN=filter(Bf,Af,NN);
                        %Concatenate the results
                        Zz=[Zz,ZZ];Ee=[Ee,EE];Nn=[Nn,NN];
                        NameZ=[NameZ,{list(3).name(1:4)}];  
                    end
                end
                %Selection of the station name and coordinates
                XS=[];YS=[];ZS=[];%Initialize the coordinates vectors
                %Loop through the number of station used
                for bb=1:length(NameZ)
                    %Selection of the station name and coordinates
                    idx=strcmp(station,NameZ{bb});
                    xss=xs(idx);
                    yss=ys(idx);
                    zss=zs(idx);
                    XS=[XS;xss];
                    YS=[YS;yss];
                    ZS=[ZS;zss];
                end
                %Picking mode 
                w1=round(w_sem/dt);ws=round(w_sem/2/dt);
                if delim==3
                    %Manual discrete mode
                    %Plot the seismic trace of reference 
                    new_figure=figure('units','normalized','outerposition',[0 0 1 1],'Name','Select window analysis','NumberTitle','off','MenuBar','none');
                    YY=Zz(:,staidx);
                    YY=YY./max(YY);
                    plot(YY)
                    ylim([-2 2])
                    xlim([1 length(YY)])
                    title([NameZ{staidx} '-' name])
                    grid on
                    grid minor
                    %Select manually the picking time
                    [L]=ginput(1);
                    trig=round(L(1,1));
                    if length(trig)>1;trig=trig(1);end
                    close(new_figure)
                elseif delim==2
                    %Automatic discrete mode
                    YY=Zz(:,staidx);
                    YY=YY./max(YY);
                    %Search the down onset 
                    [~,posmax]=max(YY);
                    sgn=YY(posmax-w1+1:posmax);
                    valmin=min(sgn);
                    trig=find(valmin==YY);
                    if length(trig)>1;trig=trig(1);end
                elseif delim==1
                    %Automatic continous mode
                    k=2;%Exclude the analysis windows on the edges of the trace
                    i1=0;i2=0;nt=size(Zz,1);trig=[];
                    %Move the analysis window through the seismic trace of reference
                    while i2 < nt-2*w1 || i2== nt-2*w1
                        k=k+1;
                        i1=(1+(k-1)*w1); i2=i1+w1-1;
                        %Select the end-limit (i2) as picking time
                        trig=[trig;i2];
                    end
                end
                %Concatenate the results
                TRIG=[TRIG;{trig}];
                XSS=[XSS;{XS}];
                YSS=[YSS;{YS}];
                ZSS=[ZSS;{ZS}];
                Z=[Z;{Zz}];
                E=[E;{Ee}];
                N=[N;{Nn}];
                Namemain=[Namemain;name];
            end
        end
    catch
        %Error assesment
        return
    end
end
%==========================================================================
%%%%Plot of XY section of Radial Semblance in UTM system %%%%%%%%%%%%%%%%%%
function ploting1(ii,ax)
    %Acquisition of the Grid limits
    xmin=str2num(get(h.xmin,'string'));
    xmax=str2num(get(h.xmax,'string'));
    ymin=str2num(get(h.ymin,'string'));
    ymax=str2num(get(h.ymax,'string'));
    zmin=str2num(get(h.zmin,'string'));
    zmax=str2num(get(h.zmax,'string'));
    loc=cell2mat(GRID(ii));%Select the ii-th 3D Radial Semblance grid
    %Select the ii-th maximum values of analysis grid (x,y,z, Radial Semblance value)
    ix1=LOCAL(ii,1);
    ix2=LOCAL(ii,2);
    ix3=LOCAL(ii,3);
    valor=LOCAL(ii,4);
    %Select the ii-th filename
    name=Name(ii,:);
    % Set the title with the analysis results
    tit1=cellstr(name);
    tit2=cellstr(sprintf('Max point (%.3f) : Nx=%.2f, Ny=%.2f, Nz=%.2f \n',valor,xq(ix2),yq(ix1),zq(ix3)));
    tit3=cellstr(sprintf('Error Max point  : Ex=%.3f, Ey=%.3f, Ez=%.3f \n',ERROR(ii,1),ERROR(ii,2),ERROR(ii,3)));
    tit=[tit1;tit2;tit3];
    cla(ax)%Clear XY plane axes
    %Plot of XY section of Radial Semblance  
    C=squeeze(loc(:,:,ix3)); 
    hold(ax,'on') 
    surf(ax,squeeze(Xq(:,:,ix3)),squeeze(Yq(:,:,ix3)),squeeze(Zq(:,:,ix3)),C);
    shading(ax,'interp')
    caxis(ax,[0 valor])
    az=0;el=90;
    ylabel(ax,'Y (km)')
    xlabel(ax,'X (km)')
    view(ax,az,el)
    colorbar(ax)
    colormap(ax,'parula');
    %Plotting the stations used
    scatter3(ax,xs,ys,zs,'ok','filled')
    xlim(ax,[xmin xmax]);ylim(ax,[ymin ymax])
    title(ax,tit)
    hold(ax,'off')
end

%==========================================================================
%%%%Plot of XY section of Radial Semblance in wgs84 system %%%%%%%%%%%%%%%%
function ploting1_deg(ii,ax)
    %Acquisition of the Grid limits
    xmin=str2num(get(h.xmin,'string'));
    xmax=str2num(get(h.xmax,'string'));
    ymin=str2num(get(h.ymin,'string'));
    ymax=str2num(get(h.ymax,'string'));
    min=str2num(get(h.zmin,'string'));
    zmax=str2num(get(h.zmax,'string'));
    %Convert UTM coordinates (Latitude, Longitude) into WGS84 coordinates 
    [ymin,xmin]=utm2deg(xmin.*1000,ymin.*1000,utm_zone(1,:));
    [ymax,xmax]=utm2deg(xmax.*1000,ymax.*1000,utm_zone(1,:));
    loc=cell2mat(GRID(ii));%Select the ii-th 3D Semblance grid
    %Select the ii-th maximum values of analysis grid (x,y,z,Radial Semblance value)
    ix1=LOCAL(ii,1);
    ix2=LOCAL(ii,2);
    ix3=LOCAL(ii,3);
    valor=LOCAL(ii,4);
    %Select the ii-th filename
    name=Name(ii,:);
    %Conversion factors in analysis error assessment
    ferry=yyq(2)-yyq(1);ferrx=xxq(2)-xxq(1);
    % Set the title with the analysis results
    tit1=cellstr(name);
    tit2=cellstr(sprintf('Max point (%.3f) : Lat=%.3f, Lon=%.3f, Elv=%.3f \n',valor,yyq(ix1),xxq(ix2),zzq(ix3)));
    tit3=cellstr(sprintf('Error Max point  : Ex=%.4f, Ey=%.4f, Ez=%.3f \n',(ERROR(ii,2)*ferry)/100,(ERROR(ii,1)*ferrx)/100,ERROR(ii,3)));
    tit=[tit1;tit2;tit3];
    cla(ax)%Clear XY plane axes
    %Plot of XY section of Radial Semblance 
    C=squeeze(loc(:,:,ix3)); 
    hold(ax,'on') 
    surf(ax,squeeze(XXq(:,:,ix3)),squeeze(YYq(:,:,ix3)),squeeze(ZZq(:,:,ix3)),C);
    shading(ax,'interp')
    caxis(ax,[0 valor])
    az=0;el=90;
    ylabel(ax,'Latitude (decimal degrees)')
    xlabel(ax,'Longitude (decimal degrees)')
    view(ax,az,el)
    colorbar(ax)
    colormap(ax,'parula');
    %Plotting the stations used
    scatter3(ax,LON,LAT,ELE./1000,'ok','filled')
    xlim(ax,[xmin xmax]);ylim(ax,[ymin ymax])
    title(ax,tit)
    hold(ax,'off')
end

%==========================================================================
%%%%Plot of YZ section of Radial Semblance in UTM system %%%%%%%%%%%%%%%%%%
function ploting2(ii,ax)
    %Acquisition of the Grid limits
    xmin=str2num(get(h.xmin,'string'));
    xmax=str2num(get(h.xmax,'string'));
    ymin=str2num(get(h.ymin,'string'));
    ymax=str2num(get(h.ymax,'string'));
    zmin=str2num(get(h.zmin,'string'));
    zmax=str2num(get(h.zmax,'string'));
    %Select the ii-th 3D Radial Semblance grid
    loc=cell2mat(GRID(ii));
    %Select the ii-th maximum values of analysis grid (x,y,z,Radial Semblance value)
    ix2=LOCAL(ii,2);
    valor=LOCAL(ii,4);
    cla(ax)%Clear YZ plane axes
    %Plot of YZ section of Radial Semblance 
    colormap(ax,'parula');
    C=squeeze(loc(:,ix2,:)); 
    surf(ax,squeeze(Xq(:,ix2,:)),squeeze(Yq(:,ix2,:)),squeeze(Zq(:,ix2,:)),C);
    shading(ax,'interp')
    caxis(ax,[0 valor])
    az=90;el=0;
    ylabel(ax,'Y (km)')
    zlabel(ax,'Z (km)')
    view(ax,az,el)
    colorbar(ax)
    zlim(ax,[zmin zmax]);ylim(ax,[ymin ymax])
end

%==========================================================================
%%%%Plot of YZ section of Radial Semblance in wgs84 system %%%%%%%%%%%%%%%%
function ploting2_deg(ii,ax)
    %Acquisition of the Grid limits
    xmin=str2num(get(h.xmin,'string'));
    xmax=str2num(get(h.xmax,'string'));
    ymin=str2num(get(h.ymin,'string'));
    ymax=str2num(get(h.ymax,'string'));
    zmin=str2num(get(h.zmin,'string'));
    zmax=str2num(get(h.zmax,'string'));
    %Convert UTM coordinates (Latitude, Longitude) into WGS84 coordinates 
    [ymin,xmin]=utm2deg(xmin.*1000,ymin.*1000,utm_zone(1,:));
    [ymax,xmax]=utm2deg(xmax.*1000,ymax.*1000,utm_zone(1,:));
    %Select the ii-th 3D Radial Semblance grid
    loc=cell2mat(GRID(ii));
    %Select the ii-th maximum values of analysis grid (x,y,z,Radial Semblance value)
    ix2=LOCAL(ii,2);
    valor=LOCAL(ii,4);
    cla(ax)%Clear YZ plane axes
    %Plot of YZ section of Radial Semblance
    colormap(ax,'parula');
    C=squeeze(loc(:,ix2,:)); 
    surf(ax,squeeze(XXq(:,ix2,:)),squeeze(YYq(:,ix2,:)),squeeze(ZZq(:,ix2,:)),C);
    shading(ax,'interp')
    caxis(ax,[0 valor])
    az=90;el=0;
    ylabel(ax,'Latitude (decimal degrees)')
    zlabel(ax,'Elevation(km)')
    view(ax,az,el)
    colorbar(ax)
    zlim(ax,[zmin zmax]);ylim(ax,[ymin ymax])
end

%==========================================================================
%%%%Plot of XZ section of Radial Semblance in UTM system %%%%%%%%%%%%%%%%%%%%%%%%%
function ploting3(ii,ax)
    %Acquisition of the Grid limits
    xmin=str2num(get(h.xmin,'string'));
    xmax=str2num(get(h.xmax,'string'));
    ymin=str2num(get(h.ymin,'string'));
    ymax=str2num(get(h.ymax,'string'));
    zmin=str2num(get(h.zmin,'string'));
    zmax=str2num(get(h.zmax,'string'));
    %Select the ii-th 3D Radial Semblance grid
    loc=cell2mat(GRID(ii));
    %Select the ii-th maximum values of analysis grid (x,y,z,Radial emblance value)
    ix1=LOCAL(ii,1);
    valor=LOCAL(ii,4);
    cla(ax)%Clear XZ plane axes
    %Plot of XZ section of Radial Semblance
    colormap(ax,'parula');
    C=squeeze(loc(ix1,:,:)); 
    surf(ax,squeeze(Xq(ix1,:,:)),squeeze(Yq(ix1,:,:)),squeeze(Zq(ix1,:,:)),C);
    shading(ax,'interp')
    caxis(ax,[0 valor])
    az=0;el=0;
    xlabel(ax,'X (km)')
    zlabel(ax,'Z (km)')
    view(ax,az,el)
    colorbar(ax)
    xlim(ax,[xmin xmax]);zlim(ax,[zmin zmax])
end

%==========================================================================
%%%%Plot of XZ section of Radial Semblance in wgs84 system %%%%%%%%%%%%%%%%
function ploting3_deg(ii,ax)
    %Acquisition of the Grid limits
    xmin=str2num(get(h.xmin,'string'));
    xmax=str2num(get(h.xmax,'string'));
    ymin=str2num(get(h.ymin,'string'));
    ymax=str2num(get(h.ymax,'string'));
    zmin=str2num(get(h.zmin,'string'));
    zmax=str2num(get(h.zmax,'string'));
    %Convert UTM coordinates (Latitude, Longitude) into WGS84 coordinates 
    [ymin,xmin]=utm2deg(xmin.*1000,ymin.*1000,utm_zone(1,:));
    [ymax,xmax]=utm2deg(xmax.*1000,ymax.*1000,utm_zone(1,:));
    %Select the ii-th 3D Radial Semblance grid
    loc=cell2mat(GRID(ii));
    %Select the ii-th maximum values of analysis grid (x,y,z,Radial Semblance value)
    ix1=LOCAL(ii,1);
    valor=LOCAL(ii,4);
    cla(ax)%Clear XZ plane axes
    %Plot of XZ section of Radial Semblance
    colormap(ax,'parula');
    C=squeeze(loc(ix1,:,:)); 
    surf(ax,squeeze(XXq(ix1,:,:)),squeeze(YYq(ix1,:,:)),squeeze(ZZq(ix1,:,:)),C);
    shading(ax,'interp')
    caxis(ax,[0 valor])
    az=0;el=0;
    xlabel(ax,'Longitude (decimal degrees)')
    zlabel(ax,'Elevation (m)')
    view(ax,az,el)
    colorbar(ax)
    xlim(ax,[xmin xmax]);zlim(ax,[zmin zmax])
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
%%%%Refresh the XYZ sections axes through the pre button%%%%%%%%%%%%%%%%%%%
function loca1(~,~)
   %Error control
   if isempty(LOCAL);set(text1, 'String','Error. Localization matrix is empty.');drawnow;return;end
   %Decrease the slidval value
   slidval=slidval-1;if slidval<=0;slidval=1;end
   %Plot the new blocks on the basis of the type of coordinates system (UTM/wgs84 system) 
   switch params.coord_sys
       case 'm'
           %Plot the new block on the XY section
           ploting1(slidval,axes1);
           %Plot the new block on the YZ section
           ploting2(slidval,axes2);
           %Plot the new block on the XZ section
           ploting3(slidval,axes3);
       case 'degrees'
           %Plot the new block on the XY section
           ploting1_deg(slidval,axes1);
           %Plot the new block on the YZ section
           ploting2_deg(slidval,axes2);
           %Plot the new block on the XZ section
           ploting3_deg(slidval,axes3);
   end        
end

%==========================================================================
%%%%Refresh the XYZ sections axes through the next button%%%%%%%%%%%%%%%%%%
function loca2(~,~)
   len=size(ff,2);
   %Error control
   if isempty(LOCAL);set(text1, 'String','Error. Localization matrix is empty.');return;end
   %Increase the slidval value
   slidval=slidval+1;if slidval>len;slidval=len;end
   %Plot the new blocks on the basis of the type of coordinates system (UTM/wgs84 system) 
   switch params.coord_sys
       case 'm'
           %Plot the new block on the XY section
           ploting1(slidval,axes1);
           %Plot the new block on the YZ section
           ploting2(slidval,axes2);
           %Plot the new block on the XZ section
           ploting3(slidval,axes3);
       case 'degrees'
           %Plot the new block on the XY section
           ploting1_deg(slidval,axes1);
           %Plot the new block on the YZ section
           ploting2_deg(slidval,axes2);
           %Plot the new block on the XZ section
           ploting3_deg(slidval,axes3);
   end         
end
%==========================================================================
%%%%Set the limits of the x-axis on the XYZ axes%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setlimits1(~,~)
    try
        %Manual settings
        %Acquisition limits
        x1_lim=str2num(get(h.x1_lim,'string'));
        x2_lim=str2num(get(h.x2_lim,'string'));  
        xlim(axes1,[x1_lim x2_lim])
        xlim(axes2,[x1_lim x2_lim])
        xlim(axes3,[x1_lim x2_lim])
    catch
        %Default settings
        if ~isempty(LOCAL)
            switch params.coord_sys
                case 'm'
                    %Reset limits (UTM system)
                    xlim(axes1,[min(xq) max(xq)])
                    xlim(axes2,[min(xq) max(xq)])
                    xlim(axes3,[min(xq) max(xq)])
                case 'degrees'
                    %Reset limits (wgs84)
                    xlim(axes1,[min(xxq) max(xxq)])
                    xlim(axes2,[min(xxq) max(xxq)])
                    xlim(axes3,[min(xxq) max(xxq)])
            end
        end
    end
end

%==========================================================================
%%%%Set the limits of the y-axis on the XYZ axes%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setlimits2(~,~)
    try
        %Manual settings
        %Acquisition limits
        y1_lim=str2num(get(h.y1_lim,'string'));
        y2_lim=str2num(get(h.y2_lim,'string'));
        ylim(axes1,[y1_lim y2_lim])
        ylim(axes2,[y1_lim y2_lim])
        ylim(axes3,[y1_lim y2_lim])
    catch
        %Default settings
        if ~isempty(LOCAL)
            switch params.coord_sys
                case 'm'
                    %Reset limits (UTM system)
                    ylim(axes1,[min(yq) max(yq)])
                    ylim(axes2,[min(yq) max(yq)])
                    ylim(axes3,[min(yq) max(yq)])
                case 'degrees'
                    %Reset limits (wgs84)
                    ylim(axes1,[min(yyq) max(yyq)])
                    ylim(axes2,[min(yyq) max(yyq)])
                    ylim(axes3,[min(yyq) max(yyq)])
            end
        end
    end
end

%==========================================================================
%%%%Set the limits of the z-axis on the XYZ axes%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setlimits3(~,~)
    try
        %Manual settings
        %Acquisition limits
        z1_lim=str2num(get(h.z1_lim,'string'));
        z2_lim=str2num(get(h.z2_lim,'string'));
        zlim(axes1,[z1_lim z2_lim])
        zlim(axes2,[z1_lim z2_lim])
        zlim(axes3,[z1_lim z2_lim])
    catch
        %Default settings
        if ~isempty(LOCAL)
            switch params.coord_sys
                case 'm'
                    %Reset limits (UTM system)
                    zlim(axes1,[min(zq) max(zq)])
                    zlim(axes2,[min(zq) max(zq)])
                    zlim(axes3,[min(zq) max(zq)])
                case 'degrees'
                    %Reset limits (wgs84)
                    zlim(axes1,[min(zzq) max(zzq)])
                    zlim(axes2,[min(zzq) max(zzq)])
                    zlim(axes3,[min(zzq) max(zzq)])
            end
        end
    end
end

%==========================================================================
%%%%Control and reset automatically analysis parameters%%%%%%%%%%%%%%%%%%%%
function deflimits(~,~)
        %Acquisition analysis parameters
        w_sem=str2num(get(h.w_sem,'string'));
        vel_sem=str2num(get(h.vel_sem,'string'));
        freq_sem=get(h.freq_sem,'string');
        xmin=str2num(get(h.xmin,'string'));
        xmax=str2num(get(h.xmax,'string'));
        ymin=str2num(get(h.ymin,'string'));
        ymax=str2num(get(h.ymax,'string'));
        zmin=str2num(get(h.zmin,'string'));
        zmax=str2num(get(h.zmax,'string'));
        
        %Control and eventually reset parameters
        if (size(w_sem,1)==0);set(h.w_sem,'string',params.w_sem);end
        if (size(vel_sem,1)==0);set(h.vel_sem,'string',params.vel_sem);end
        if (size(freq_sem,2)==0);set(h.freq_sem,'string',params.freq_sem);end
        if (size(xmin,1)==0);set(h.xmin,'string',params.xmin);end
        if (size(ymin,1)==0);set(h.ymin,'string',params.ymin);end
        if (size(zmin,1)==0);set(h.zmin,'string',params.zmin);end
        if (size(xmax,1)==0);set(h.xmax,'string',params.xmax);end
        if (size(ymax,1)==0);set(h.ymax,'string',params.ymax);end
        if (size(zmax,1)==0);set(h.zmax,'string',params.zmax);end
end

%==========================================================================
%%%%%%%%%%%%%%%%%%%%% Save the main plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save1(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Acquisition and validation of the Grid limits
    xmin=str2num(get(h.xmin,'string'));
    xmax=str2num(get(h.xmax,'string'));
    ymin=str2num(get(h.ymin,'string'));
    ymax=str2num(get(h.ymax,'string'));
    zmin=str2num(get(h.zmin,'string'));
    zmax=str2num(get(h.zmax,'string'));
    if xmin<=0;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if xmax<=0 | xmax<=xmin;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if ymin<=0;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if ymax<=0 | ymax<=ymin;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if zmin<=0;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if zmax<=0 | zmax<=zmin;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if params.step_sem<=0;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    try
        %Set the block of results
        ii=slidval;
        %Select the ii-th 3D Radial Semblance grid
        loc=cell2mat(GRID(ii));
        %Select the ii-th maximum values of analysis grid (x,y,z,Radial Semblance value)
        ix1=LOCAL(ii,1);
        ix2=LOCAL(ii,2);
        ix3=LOCAL(ii,3);
        valor=LOCAL(ii,4);
        %Select the ii-th filename
        name=Name(ii,:);
        %Plot the results on the basis of the type of coordinates system (UTM/wgs84 system) 
        switch params.coord_sys
            case 'm'
                new_figure=figure('units','normalized','outerposition',[0 0 1 1],'Name','Radial Semblance Grid','NumberTitle','off');
                subplot(2,2,1)
                %Plot of XY section of Radial Semblance
                C=squeeze(loc(:,:,ix3));
                hold('on') 
                surf(squeeze(Xq(:,:,ix3)),squeeze(Yq(:,:,ix3)),squeeze(Zq(:,:,ix3)),C);
                shading('interp')
                caxis([0 valor])
                az=0;el=90;
                ylabel('Y (km)')
                xlabel('X (km)')
                view(az,el)
                colorbar
                colormap('parula');
                grid on 
                grid minor
                %Plotting the stations used
                scatter3(xs,ys,zs,'ok','filled')
                xlim([xmin xmax]);ylim([ymin ymax])
                title(tit)
                hold('off')
                subplot(2,2,2)
                %Plot of YZ section of Radial Semblance
                colormap('parula');
                C=squeeze(loc(:,ix2,:)); 
                surf(squeeze(Xq(:,ix2,:)),squeeze(Yq(:,ix2,:)),squeeze(Zq(:,ix2,:)),C);
                shading('interp')
                caxis([0 valor])
                az=90;el=0;
                ylabel('Y (km)')
                zlabel('Z (km)')
                view(az,el)
                colorbar
                grid on 
                grid minor
                zlim([zmin zmax]);ylim([ymin ymax])
                subplot(2,2,3)
                %Plot of XZ section of Radial Semblance
                colormap('parula');
                C=squeeze(loc(ix1,:,:)); 
                surf(squeeze(Xq(ix1,:,:)),squeeze(Yq(ix1,:,:)),squeeze(Zq(ix1,:,:)),C);
                shading('interp')
                caxis([0 valor])
                az=0;el=0;
                xlabel('X (km)')
                zlabel('Z (km)')
                view(az,el)
                colorbar
                grid on 
                grid minor
                xlim([xmin xmax]);zlim([zmin zmax])
            case 'degrees'
                %Convert UTM coordinates (Latitude, Longitude) into WGS84 coordinates 
                [ymin,xmin]=utm2deg(xmin.*1000,ymin.*1000,utm_zone(1,:));
                [ymax,xmax]=utm2deg(xmax.*1000,ymax.*1000,utm_zone(1,:));
                new_figure=figure('units','normalized','outerposition',[0 0 1 1],'Name','Radial Semblance Grid','NumberTitle','off');
                subplot(2,2,1)
                %Plot of XY section of Radial Semblance
                C=squeeze(loc(:,:,ix3));
                hold('on') 
                surf(squeeze(XXq(:,:,ix3)),squeeze(YYq(:,:,ix3)),squeeze(ZZq(:,:,ix3)),C);
                shading('interp')
                caxis([0 valor])
                az=0;el=90;
                ylabel('Latitude (decimal degrees)')
                xlabel('Longitude (decimal degrees)')
                view(az,el)
                colorbar
                colormap('parula');
                grid on 
                grid minor
                %Plotting the stations used
                scatter3(LON,LAT,ELE*1000,'ok','filled')
                xlim([xmin xmax]);ylim([ymin ymax])
                title(tit)
                hold('off')
                subplot(2,2,2)
                %Plot of YZ section of Radial Semblance
                colormap('parula');
                C=squeeze(loc(:,ix2,:)); 
                surf(squeeze(XXq(:,ix2,:)),squeeze(YYq(:,ix2,:)),squeeze(ZZq(:,ix2,:)),C);
                shading('interp')
                caxis([0 valor])
                az=90;el=0;
                ylabel('Latitude (decimal degrees)')
                zlabel('Elevation (km)')
                view(az,el)
                colorbar
                grid on 
                grid minor
                zlim([zmin zmax]);ylim([ymin ymax])
                subplot(2,2,3)
                %Plot of XZ section of Radial Semblance
                colormap('parula');
                C=squeeze(loc(ix1,:,:)); 
                surf(squeeze(XXq(ix1,:,:)),squeeze(YYq(ix1,:,:)),squeeze(ZZq(ix1,:,:)),C);
                shading('interp')
                caxis([0 valor])
                az=0;el=0;
                xlabel('Longitude (decimal degrees)')
                zlabel('Elevation (km)')
                view(az,el)
                colorbar
                grid on 
                grid minor
                xlim([xmin xmax]);zlim([zmin zmax])          
        end
    catch
        %Error assessment
        if isempty(LOCAL);set(text1, 'String','Error. Localization matrix is empty.');drawnow;return;end
        if size(loc,1)~=size(Xq,1)|size(loc,2)~=size(Xq,2)|size(loc,3)~=size(Xq,3);set(text1, 'String','Error. Incoherent dimensions of the Localization array.');drawnow;return;end
    end
end

%==========================================================================
%%%%%%%%%%%%%%%%%%%%% Plot the maximum probability region%%%%%%%%%%%%%%%%%%
function save2(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Acquisition and validation of the Grid limits
    xmin=str2num(get(h.xmin,'string'));
    xmax=str2num(get(h.xmax,'string'));
    ymin=str2num(get(h.ymin,'string'));
    ymax=str2num(get(h.ymax,'string'));
    zmin=str2num(get(h.zmin,'string'));
    zmax=str2num(get(h.zmax,'string'));
    if xmin<=0;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if xmax<=0 | xmax<=xmin;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if ymin<=0;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if ymax<=0 | ymax<=ymin;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if zmin<=0;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if zmax<=0 | zmax<=zmin;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if params.step_sem<=0;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    try
        %Set the block of results
        ii=slidval;
        %Select the ii-th 3D Radial Semblance grid
        loc=GRID{ii};
        %Replace the NaN values
        idx=isnan(loc);
        loc(idx)=0;
        %Select the ii-th maximum values of analysis grid (x,y,z,Radial Semblance value)
        ix1=LOCAL(ii,1);
        ix2=LOCAL(ii,2);
        ix3=LOCAL(ii,3);
        valor=LOCAL(ii,4);
        %Select the ii-th filename
        name=Name(ii,:);
        S=loc;
        %Plot the results on the basis of the type of coordinates system (UTM/wgs84 system) 
        switch params.coord_sys
            case 'm'
                %Construction of  the Search grid (no topographic correction) and the DEM grid (UTM system)
                [XXXq,YYYq,ZZZq,Nx,Ny,Nz,xxxq,yyyq,zzzq,A,X,Y]=DEMconst('mapfile_m.mat','no',xmin,xmax,ymin,ymax,zmin,zmax,params.step_sem);
                zmin=min(min(A));zmax=max(max(A));
                %Plot the 10% of the Semblance volume around the maximum value
                hnew=slice(XXXq,YYYq,ZZZq,S,xxxq(ix2),yyyq(ix1),zzzq(ix3));
                new_figure=figure('Name','Volume plot','NumberTitle','off');
                caxis([0 valor])
                colormap white
                set(hnew,'FaceColor','None','EdgeColor','none')
                factor_isoval=9.0;
                isoval = max(S(:))/10*factor_isoval;
                hnew=patch(isosurface(XXXq,YYYq,ZZZq,S,isoval),...
                    'FaceColor','red',...
                    'FaceAlpha',0.3,...
                    'EdgeColor','Red',...
                    'AmbientStrength',.5,...
                    'SpecularStrength',.7,...
                    'DiffuseStrength',.4);
                isonormals(XXXq,YYYq,ZZZq,S,hnew)
                patch(isocaps(XXXq,YYYq,ZZZq,S,isoval),...
                    'FaceColor','interp',...
                    'EdgeColor','none')
                hold on 
                %Plotting the location of the maximum Radial Semblance value
                scatter3 (xxxq(ix2),yyyq(ix1),zzzq(ix3),'ok','filled')
                %Plotting the stations used
                scatter3( xs,ys,zs,'or','filled')
                xlabel('X (m)')  
                ylabel('Y (m)')
                zlabel('Altitude (km)')
                %Plot the DEM 
                contour3(X,Y,A,(zmin:params.step_sem:zmax),'k')
                xlim([xmin xmax])
                ylim([ymin ymax])
                zlim([zmin zmax])  
                title(tit)
            case 'degrees'
                %Construction of  the Search grid (no topographic correction) and the DEM grid (wgs84 system)
                [XXXq,YYYq,ZZZq,xxxq,yyyq,zzzq,X,Y,A]=DEMconstDeg('mapfile_deg.mat','no',xmin,xmax,ymin,ymax,zmin,zmax,params.step_sem,utm_zone);
                %Convert UTM coordinates (Latitude, Longitude) into WGS84 coordinates 
                [ymin,xmin]=utm2deg(xmin.*1000,ymin.*1000,utm_zone(1,:));
                [ymax,xmax]=utm2deg(xmax.*1000,ymax.*1000,utm_zone(1,:));
                zmin=min(min(A));zmax=max(max(A));
                %Plot the 10% of the Radial Semblance volume around the maximum value
                new_figure=figure('Name','Volume plot','NumberTitle','off');
                hnew=slice(XXXq,YYYq,ZZZq,S,xxxq(ix2),yyyq(ix1),zzzq(ix3));
                caxis([0 valor])
                colormap white
                set(hnew,'FaceColor','None','EdgeColor','none')
                factor_isoval=9.0;
                isoval = max(S(:))/10*factor_isoval;
                hnew=patch(isosurface(XXXq,YYYq,ZZZq,S,isoval),...
                    'FaceColor','red',...
                    'FaceAlpha',0.3,...
                    'EdgeColor','Red',...
                    'AmbientStrength',.5,...
                    'SpecularStrength',.7,...
                    'DiffuseStrength',.4);
                isonormals(XXXq,YYYq,ZZZq,S,hnew)
                patch(isocaps(XXXq,YYYq,ZZZq,S,isoval),...
                    'FaceColor','interp',...
                    'EdgeColor','none')
                hold on 
                %Plotting the location of the maximum Radial Semblance value
                scatter3 (xxxq(ix2),yyyq(ix1),zzzq(ix3),'ok','filled')
                %Plotting the stations used
                scatter3(LON,LAT,ELE/1000,'or','filled')
                xlabel('Longitude (decimal degrees)')  
                ylabel('Latitude (decimal degrees)')
                zlabel('Elevation (km)')
                %Plot the DEM 
                contour3(X*1000,Y*1000,A,(zmin:params.step_sem:zmax),'k')
                xlim([xmin xmax])
                ylim([ymin ymax]) 
                title(tit)
        end
    catch
        %Error assessment
        if isempty(LOCAL);set(text1, 'String','Error. Localization matrix is empty.');drawnow;return;end
        if size(loc,1)~=size(XXXq,1)|size(loc,2)~=size(XXXq,2)|size(loc,3)~=size(XXXq,3);set(text1, 'String','Error. Incoherent dimensions of the Localization array.');drawnow;return;end       
    end
end

%==========================================================================
%%%%%%%%%%%%%%%%%% Plot all localizations on the DEM%%%%%%%%%%%%%%%%%%%%%%%
function save3(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Acquisition and validation of the Grid limits
    xmin=str2num(get(h.xmin,'string'));
    xmax=str2num(get(h.xmax,'string'));
    ymin=str2num(get(h.ymin,'string'));
    ymax=str2num(get(h.ymax,'string'));
    zmin=str2num(get(h.zmin,'string'));
    zmax=str2num(get(h.zmax,'string'));
    if xmin<=0;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if xmax<=0 | xmax<=xmin;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if ymin<=0;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if ymax<=0 | ymax<=ymin;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if zmin<=0;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if zmax<=0 | zmax<=zmin;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    if params.step_sem<=0;set(text1, 'String','Error. Invalid Grid value.');drawnow;return;end
    try
        %Initialize the DEM variables
        X=[];Y=[];A=[];
        %Calculation of the source locations
        x=xq(LOCAL(:,2));
        y=yq(LOCAL(:,1));
        zz=zq(LOCAL(:,3));
        valor=LOCAL(:,4)';
        %Active/Deactive the visualization of the analysis errors
        wx=ERROR(:,1);
        wy=ERROR(:,2);
        wz=ERROR(:,3);
        if emod=='y' & ~isnan(ERROR)
            %Calculation of the location errors
            cond=wx==0;
            wx(cond)=params.step_sem;
            cond=wy==0;
            wy(cond)=params.step_sem;
            cond=wz==0;
            wz(cond)=params.step_sem;
            Ex=(sum(x.*(wx)))./sum(wx);
            Ey=(sum(y.*(wy)))./sum(wy);
            Ez=(sum(zz.*(wz)))./sum(wz);

            VEx=(sum(ERROR(:,1).*(wx)))./sum(wx);
            VEy=(sum(ERROR(:,2).*(wy)))./sum(wy);
            VEz=(sum(ERROR(:,3).*(wz)))./sum(wz);
            %Set the dimension of the points on the basis of the analysis error
            vv=sqrt(wx.^2+wy.^2+wz.^2);
        else
            %Set the standard dimension of the points on the basis of the analysis error
            vv=ones(1,length(x))*0.5;
        end
        %Plot the results on the basis of the type of coordinates system (UTM/wgs84 system) 
        switch params.coord_sys
            case 'm'
                %Load the DEM variables
                load('mapfile_m.mat','X','Y','A')
                zmin=min(min(A));zmax=max(max(A));
                new_figure=figure('Name','Scatter Plot','NumberTitle','off');
                %Plot the DEM
                ax1=axes;
                contour3(ax1,X,Y,A,(zmin:params.step_sem:zmax),'k')
                xlim(ax1,[xmin xmax])
                ylim(ax1,[ymin ymax])
                zlim(ax1,[zmin zmax])
                xlabel(ax1,'X (km)')  
                ylabel(ax1,'Y (km)')
                zlabel(ax1,'Altitude (km)')
                view(ax1,[0 90])
                ax2=axes;
                xlim(ax2,[xmin xmax])
                ylim(ax2,[ymin ymax])
                zlim(ax2,[zmin 10])
                linkaxes([ax2,ax1])
                % Plot the source locations
                scatter3 (ax2,x,y,zz,vv*200,valor','filled')
                hold(ax2,'on')
                %Plot the stations used
                scatter3(ax2, xs,ys,zs,'dk','filled')
                xlim(ax2,[xmin xmax])
                ylim(ax2,[ymin ymax])
                zlim(ax2,[zmin 10])
                ax2.Visible = 'off';
                ax2.XTick = [];
                ax2.YTick = [];
                cb=colorbar(ax2);
                cb.Position = cb.Position + [.076 0 0 0];
                colormap(ax2,'parula')
                view(ax2,[0 90])
                title(tit) 
            case 'degrees'
                %Load the DEM variables
                load('mapfile_deg.mat','X','Y','A')
                len=size(ff,2);
                %Convert UTM coordinates (Latitude, Longitude) into WGS84 coordinates 
                [ymin,xmin]=utm2deg(xmin.*1000,ymin.*1000,utm_zone(1,:));
                [ymax,xmax]=utm2deg(xmax.*1000,ymax.*1000,utm_zone(1,:));
                zmin=min(min(A));zmax=max(max(A));
                [y,x]=utm2deg(x.*1000,y.*1000,repmat(utm_zone(1,:),len,1));
                new_figure=figure('Name','Scatter Plot','NumberTitle','off');
                ax1=axes;
                %Plot the DEM 
                contour3(ax1,X,Y,A,(zmin:params.step_sem:zmax),'k')
                xlim(ax1,[xmin xmax])
                ylim(ax1,[ymin ymax])
                zlim(ax1,[zmin zmax])
                xlabel(ax1,'Longitude (decimal degrees)')  
                ylabel(ax1,'Latitude (decimal degrees)')
                zlabel(ax1,'Elevation (km)')
                view(ax1,[0 90])
                ax2=axes;
                xlim(ax2,[xmin xmax])
                ylim(ax2,[ymin ymax])
                zlim(ax2,[zmin zmax])
                linkaxes([ax2,ax1])
                %Plot the source locations
                scatter3 (ax2,x,y,zz,vv*200,valor','filled')
                hold(ax2,'on')
                %Plot the stations used
                scatter3(ax2,LON,LAT,ELE/1000,'dk','filled')
                xlim(ax2,[xmin xmax])
                ylim(ax2,[ymin ymax])
                zlim(ax2,[zmin zmax])
                ax2.Visible = 'off';
                ax2.XTick = [];
                ax2.YTick = [];
                cb=colorbar(ax2);
                cb.Position = cb.Position + [.076 0 0 0];
                colormap(ax2,'parula')
                view(ax2,[0 90])
                title(tit)         
        end
    catch
        %Error assessment
        if isempty(LOCAL);set(text1, 'String','Error. Localization matrix is empty.');drawnow;return;end
        if size(y,1)~=size(x,1);set(text1, 'String','Error. Incoherent dimensions of the Localization array.');drawnow;return;end
    end
end

%==========================================================================
%%%%Save the results of the Radial Semblance analysis%%%%%%%%%%%%%%%%%%%%%%
function save_sem(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Validation of the results
    if isempty(params.save);set(text1, 'String', sprintf('Error. No output directory.'));drawnow;return;end  
    if isempty(LOCAL);set(text1, 'String', sprintf('Error. Localization matrix is empty.'));drawnow;return;end
    if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
    %Control and set the Input folder
    if ~isempty(pathnam);pathname=pathname2{slidval};pathname=[pathname '/'];end
    %Search and set the Family labels
    Fam={'F1';'F2';'F3';'F4';'No Fam'};
    words=split(pathname,{'/','\'});
    for kk=1:length(Fam)
        Famil=Fam{kk};
        idx=strcmp(words,Famil);
        val=find(idx==1);
        val=idx(val);
        if val==1
            Family=words{idx};
        break
        end
    end
    %Loop through the number of Input files (ff)
    for kk=1:size(LOCAL,1)
        name=ff{kk};time=name(end-18:end-4); 
        %Create the Output filename
        switch params.comp_chan
            case 'Chan'
                name=name(end-18:end-4); 
            case 'Comp'
                name=name(end-18:end-4);    
        end
        %Compute the temporal variables 
         yy=time(1:4);jday=time(6:8);hh=time(10:11);mi=time(12:13);ss=time(14:15);
         Time=datenum(str2num(yy),0,str2num(jday)+str2num(hh)/24+str2num(mi)/(24*60)+str2num(ss)/(24*60*60));
         %Set the kk-th results
         X=xq(LOCAL(kk,2));Y=yq(LOCAL(kk,1));Z=zq(LOCAL(kk,3));Value=LOCAL(kk,4);Ex=ERROR(kk,1);Ey=ERROR(kk,2);Ez=ERROR(kk,3);
         %Create a table of the results
         data=table(Time,X,Y,Z,Value,Ex,Ey,Ez,GRID(kk,:));
         data.Properties.VariableNames={'Time','X', 'Y','Z','Semblace','ErrorX','ErrorY','ErrorZ','Grid'};
         data.Properties.VariableUnits = {'datenum','km','km','km',' ','km','km','km','kmxkmxkm'};
         %Create the Output directory
         pathout=[params.save '/RSemblance/' Family '/' jday '/'];if ~exist(pathout, 'dir');mkdir(pathout);end 
         output=[pathout name params.outputformat];
         %Save the results
         switch params.outputformat 
             case '.mat'
                 save (output,'data');
             case '.txt'
                 Grid=cell2mat(data.Grid);
                 Grid = reshape(Grid,[size(Grid,1),size(Grid,2)*size(Grid,3)]);
                 data2=array2table(Grid);
                 writetable(data(:,1:end-1),output)
                 output2=[pathout name '_GRID' params.outputformat];
                 writetable(data2,output2)
         end
    end  
    %Command message
    set(text1, 'String', sprintf(['RSemblance' params.outputformat '-Observed data are saved.']));
    drawnow
end
end    