%%%%%%%%%%%%% Calculation of the MUSIC parameters %%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params: user parameters
% avfig: MUSIC figure
% panel_color: colors of the main panel 
% fig_color: colors of the figure
% pan1: control right side panel 
% text1: Window message
% axes1,axes2,axes3: axes of the current figure
% az1_lim, az2_lim: limit of y-axis on azimuth axes
% ray1_lim, ray2_lim: limit of y-axis on ray parameter axes
% pow1_lim, pow2_lim: limit of y-axis on fk power axes
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
% fsam: sampling rate
% nyq: nyquist number
% z: amplitude vector of the seismic trace
% nt: number of samples
% startTime, endTime: temporal limits of the traces 
%--- MUSIC setting
% typeMethod: type of analysis (Continous/Discrete)
% freq0: the central frequency of the band of interest (Hz)
% nfreq: the number of frequencies of the band of interest (Hz)
% m: the haming value
% R: the R value
% pmax: the limit of the grid slowness (s/km) 
% pinc: the size step of the grid slowness (s/km) 
% ninicio: the time origin of the analysis window (s)  
% lwin: the analysis window for the MUSIC analysis (s) 
% advance: the advancement of the analysis window for the MUSIC analysis (%) 
% nsig: the number of source to use exclusively for the plotting 
% amp: amplitude vectors of the rms
% baz1,baz2,azstep: azimuthal limits and step for histograms 3D (degrees) 
% t1,t2,nbins/timestep: temporal limits and step for histograms 3D
% eigperc: the threshold value of the sources selection
% freq_rms: frequency interval analysis
% BAZI,bazi: backazimuth vector/matrix  (degrees)
% SLOW,Slow: ray parameter vector/matrix (s/km)  
% FKPOW,fkpow: fk power vector/matrix 
% INC,Inc: incidence vector/matrix (degrees)
% SXMAX/Sxmax, SYMAX/Symax: horizontal slowness vectors/matrixs (s/km) 
% Time, TIME_music: time vectors
% TMP, TMP_perc: the cumulative matrix/vector  

function musicArray()
%% Declaration global parameters
clc
global params
comp=params.component;component=[];
pmax=params.pmax;  
pinc=params.pinc; 
ninicio=params.ninicio; 
typeMethod='Con';
ff= cell(1);
mytrace=[];
slidval=1;nsig=0;fsam=[];
pathname=[];pathnam=[];pathname2=[];
TIME_music=[];
Slow=[];
bazi=[];
fkpow=[];
TMP_perc=[];
Sxmax=[];
Symax=[];
Inc=[];
xs=[];ys=[];zs=[];station={};utm_zone=[];
LAT=[];LON=[];ELE=[];
%% MUSIC figure
fig_color= [.94 .94 .94];
panel_color= [.97 .97 .97];

avfig= figure('Name', 'MUSIC', 'position',get(0,'screensize'), ...
    'NumberTitle','off','units','normalized', 'color', fig_color,'MenuBar','none');

%Window message to user
text1= uicontrol('Style','edit','FontName','Arial','Fontsize',10,'Fontweight','bold',...
                'String','Window Message','Max',5,'Min',0,'HorizontalAlignment','left',...
                'BackgroundColor', 'w','Enable','inactive',...
                'units','normalized','Position',[.05 .02 .72 .05]);


%% Axeses
% Fk-power axes
axes1= axes('position',[.05 .72 .72 .25]);
hold(axes1,'on')
     grid(axes1,'on')
     grid(axes1,'minor')
     ylim(axes1,[0 1])
     ylabel(axes1,'FK-POWER');
     colorbar(axes1)
     colormap(axes1,'parula')          
% BackAzimuth axes
axes2= axes('position',[.05 .42 .72 .25]);
hold(axes2,'on')
     grid(axes2,'on')
     grid(axes2,'minor')
     ylim(axes2,[0 360])
     ylabel(axes2,'Backazimuth (degrees)');
     colorbar(axes2)
     colormap(axes2,'parula')
% Rayparameter axes
axes3= axes('position',[.05 .12 .72 .25]);
hold(axes3,'on')
     grid(axes3,'on')
     grid(axes3,'minor')
     ylim(axes3,[0 5])
     ylabel(axes3,'RayParameter (s/km)');
     colorbar(axes3)
     colormap(axes3,'parula')

%% Control right side panel
pan1= uipanel(avfig,'visible','on','Position',[.81 .03 .185 .96],...
    'BackgroundColor', panel_color);

uipanel(avfig,'visible','on','Position',[.81 .28 .185 .01],...
    'BackgroundColor', fig_color);

%% Text boxes
%-------------------------------- MUSIC setting
uicontrol(pan1,'Style','text', 'String','MUSIC setting',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.55 .76 .4 .06],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','F0 (Hz)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .75 .3 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Nfreq',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .70 .3 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','H_point',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .65 .3 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','W (s)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .58 .3 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Adv (%)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .53 .3 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Nsig',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .48 .3 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Neig (%)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .43 .3 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','R',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .38 .3 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Vel (km/s)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .33 .3 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Xbin (h)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.55 .28 .3 .04],...
    'BackgroundColor',panel_color);

%-------------------------------- Axis setting
uicontrol(pan1,'Style','text', 'String','Axis setting',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.03 .78 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Pow1',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .75 .4 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Pow2',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .70 .4 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Az1 (deg)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .63 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Az2 (deg)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .58 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Ray1 (s/km)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .53 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Ray2 (s/km)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .48 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','t1',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .43 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','t2',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .38 .4 .04],...
    'BackgroundColor',panel_color);

%-------------------------------- MUSIC setting
h.freq0= uicontrol(pan1,'Style','edit', 'String',params.freq0,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .74 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.nfreq= uicontrol(pan1,'Style','edit', 'String',params.nfreq,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .69 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.m= uicontrol(pan1,'Style','edit', 'String',params.m,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .64 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.lwin= uicontrol(pan1,'Style','edit', 'String',params.lwin,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .59 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.advance= uicontrol(pan1,'Style','edit', 'String',params.advance,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .54 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.nsig= uicontrol(pan1,'Style','edit', 'String',params.nsig,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .49 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.eigperc= uicontrol(pan1,'Style','edit', 'String','95',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .44 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.R= uicontrol(pan1,'Style','edit', 'String',params.R,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .39 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.vel= uicontrol(pan1,'Style','edit', 'String','1.6',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .34 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);
h.nbins= uicontrol(pan1,'Style','edit', 'String','3',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.76 .29 .2 .04],...
    'BackgroundColor','w','callback',@deflimits);

%-------------------------------- Axis setting
h.pow1_lim= uicontrol(pan1,'Style','edit', 'String','0',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .74 .2 .04],...
    'BackgroundColor','w','callback',@setlimits1);
h.pow2_lim= uicontrol(pan1,'Style','edit', 'String','1',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .69 .2 .04],...
    'BackgroundColor','w','callback',@setlimits1);
h.az1_lim= uicontrol(pan1,'Style','edit', 'String','0',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .64 .2 .04],...
    'BackgroundColor','w','callback',@setlimits2);
h.az2_lim= uicontrol(pan1,'Style','edit', 'String','360',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .59 .2 .04],...
    'BackgroundColor','w','callback',@setlimits2);
h.ray1_lim= uicontrol(pan1,'Style','edit', 'String','0',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .54 .2 .04],...
    'BackgroundColor','w','callback',@setlimits3);
h.ray2_lim= uicontrol(pan1,'Style','edit', 'String','5',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .49 .2 .04],...
    'BackgroundColor','w','callback',@setlimits3);
h.t1_lim= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.18 .44 .3 .04],...
    'BackgroundColor','w','callback',@setlimits4);
h.t2_lim= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.18 .39 .3 .04],...
    'BackgroundColor','w','callback',@setlimits4);

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
                'Value',0,'Position',[.815 .37 .05 .03],'callback',@logscale);



% Editable text
uicontrol(pan1,'Style','text', 'String','MUSIC',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.03 .17 .4 .07],...
    'BackgroundColor',panel_color);

% Calculation and save button
uicontrol('style','pushbutton', 'string','Calculation','units','normalized',...
    'position',[.815 .18 .08 .0494], 'callback',@calc); 
uicontrol('style','pushbutton', 'string','Save Main Plot','units','normalized',...
    'position',[.815 .115 .08 .0494], 'callback',@save1); 
uicontrol('style','pushbutton', 'string','Save Data','units','normalized',...
    'position',[.815 .05 .08 .0494], 'callback',@save_music);
uicontrol('style','pushbutton', 'string','Ray vs. Az','units','normalized',...
    'position',[.905 .18 .08 .0494], 'callback',@save2);
uicontrol('style','pushbutton', 'string','Slow Grid','units','normalized',...
    'position',[.905 .115 .08 .0494], 'callback',@save3);
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
    TIME_music=[];Slow=[];bazi=[];fkpow=[];TMP_perc=[];Sxmax=[];Symax=[];Inc=[];%Initialize MUSIC attributes vectors
    xs=[];ys=[];zs=[];station={};LAT=[];LON=[];ELE=[];%Initialize coordinates vectors
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);cla(axes2);cla(axes3)%Clear the axes
    tt= 0;pathnam=[];pathname=[];
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
    TIME_music=[];Slow=[];bazi=[];fkpow=[];TMP_perc=[];Sxmax=[];Symax=[];Inc=[];%Initialize MUSIC attributes vectors
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
%%%%Calculation and plots the MUSIC attributes%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calc(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);cla(axes2);cla(axes3);%Clear the axes
    try
        %Acquisition and validation of the analysis parameters
        nbins=str2num(get(h.nbins,'string'));if nbins<=0;set(text1, 'String','Error. Invalid temporal bin value.');drawnow;return;end 
        freq0=str2num(get(h.freq0,'string'));if freq0<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        nfreq=str2num(get(h.nfreq,'string'));if nfreq<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        m=str2num(get(h.m,'string'));if (m~=0) & (m~=1);set(text1, 'String','Error. Invalid Hamming value.');drawnow;return;end
        lwind=str2num(get(h.lwin,'string'));if lwind<=0;set(text1, 'String','Error. Invalid analysis window.');drawnow;return;end  
        advance=str2num(get(h.advance,'string'));advance=advance/100;
        NSig1=str2num(get(h.nsig,'string'));NSig2=0;if NSig1<=0;set(text1, 'String','Error. Invalid number of sources.');drawnow;return;end
        eigperc=str2num(get(h.eigperc,'string'))./100;if eigperc<=0 | eigperc>1;set(text1, 'String','Error. Invalid eigenvalue threshold.');drawnow;return;end
        R=str2num(get(h.R,'string'));if R<=0;set(text1, 'String','Error. Invalid R value.');drawnow;return;end
        vel=str2num(get(h.vel,'string'));if vel<=0;set(text1, 'String','Error. Invalid velocity value.');drawnow;return;end
        if params.pinc<=0 | params.pmax<=0 | params.pinc>=params.pmax;set(text1, 'String','Error. Invalid slowness grid.');drawnow;return;end 
        nsig=NSig1;
        %Command window
        textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
        set(text1, 'String', textLabel); 
        drawnow
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
        nstation=length(station);%Number of station used
        %Error control
        if NSig1>length(station);set(text1, 'String', 'Error. Invalid number of sources.');drawnow;return;end
        %Reading .mat files
        [~, len]=size (ff);%Number of files selected
        TIME_music=cell(len,1);%Initialize the time matrix
        Slow=cell(len,1);%Initialize the slowness matrix
        bazi=cell(len,1);%Initialize the backazimuth matrix
        fkpow=cell(len,1);%Initialize the fk-power matrix
        TMP_perc=cell(len,1);%Initialize the cumulative matrix
        Sxmax=cell(len,1);Symax=cell(len,1);%Initialize the horizontal slowness matrix
        Inc=cell(len,1);%Initialize the incidence matrix
        x=[xs'; ys'];%Coordinates matrix
        nsta=length(station);%Number of station used
        x(1,1:nsta)=x(1,1:nsta)-mean(x(1,1:nsta));
        %%Calculate the station coordinates realtive to array center
        x(2,1:nsta)=x(2,1:nsta)-mean(x(2,1:nsta));
        %%Set weights of frequency-smoothing winodw
        %Suggested value of m is 0 (i.e. no smoothing) if CSS method is chosen
        [a,~]=haming(m); 
        %Form the focusing matrix for current frequency (see Hung and Kaveh, 1986)
        tt=eye(nsta,nsta);
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
                    %Concatenate the seismic traces through the columns
                    Z=[Z z];
                end
                dat=Z';fsam=fs;lwin=round(lwind*fsam);ninici=round(ninicio*fsam);
                %Error control
                if length(z)-lwin<=ninici | length(z)<=lwin;set(text1, 'String','Error. Invalid T0 value or window analysis.');drawnow;return;end
                frq = freq0 + (nfreq-1)*(1./(lwin*dt));
                %Error control
                if frq>=round(fsam/2);set(text1, 'String','Invalid analysis frequencies.');drawnow;return;end
                %Apply MUSIC algorithm
                [Time,SLOW,BAZI,FKPOWER,TMP,SXMAX,SYMAX,freq_rms]=musicLocalization(dat,x,ninici,fsam,dt,lwin,advance,R,nfreq,freq0,m,pmax,pinc,NSig1,NSig2,nsta,a,tt,eigperc);
                time=datetime(startTime,'ConvertFrom','datenum');%Covert in datetime format
                newtime=datenum(time+seconds(Time));%Covert in datenum format
                %Activate/Deactivate discrete analysis
                if typeMethod=='Dis'
                    %Selection of discrete values
                    imax=RMS(dat,ninici,dt,lwin,advance,freq0,m,freq_rms);
                    newtime=newtime(imax);SLOW=SLOW(imax,:);BAZI=BAZI(imax,:);FKPOWER=FKPOWER(imax,:);SXMAX=SXMAX(imax,:);SYMAX=SYMAX(imax,:);TMP=TMP(imax,:);
                end
                %Estimation of incidence
                angl=(SLOW.*vel);
                L=find(angl>1);
                angl(L)=NaN;
                inc=asind(angl);
                %Concatenate the results
                TIME_music(i)={newtime};
                Slow(i)={SLOW};
                bazi(i)={BAZI};
                fkpow(i)={FKPOWER};
                TMP_perc(i)={TMP};
                Sxmax(i)={SXMAX};
                Symax(i)={SYMAX};
                Inc(i)={inc};
            else
                %Warning message
                set(text1, 'String', 'Warning. One or more stations are no available.');
                drawnow
            end
        end
        %Error control
        if isempty(cell2mat(bazi));set(text1, 'String', 'Error. Reset Coordinates parameters.');drawnow;return;end
        %Set the block of results
        slidval=NSig1;
        %Plot the results in the fk-power axes
        ploting1(slidval,axes1)
        %Plot the results in the backazimuth axes
        ploting2(slidval,axes2)
        %Plot the results in the rayparameter axes
        ploting3(slidval,axes3)
        %Command window
        set(text1, 'String','Calculation finished.');
        drawnow
    catch
        %Error assessment
        if isempty(params.coord) | isempty(station);set(text1, 'String','Error. No coordinate file directory or No correct coordinate file.');drawnow;return;end
        if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
        if isempty(list);set(text1, 'String','Error. No files are found.');drawnow;return;end
        if isempty(z);set(text1, 'String','Error. No correct signals data.');drawnow;return;end
    end   
end
%==========================================================================
%%%%Calculation of the rms index %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imax=RMS(dat,ninici,dt,lwin,advance,freq0,m,freq_rms)
    try
        %Initialize temporal analysis variables
        iday=0;
        ihr=0;
        imin=0;
        ssec=ninici*dt;
        time_step=lwin*advance/fsam;ndat=length(dat);
        nwin=round((ndat)/(time_step*fsam));
        num_time=nwin;
        tprime=lwin;delta_time=1/fsam;
        %Set the band-pass filter coefficient
        f1=freq_rms(1);f2=freq_rms(end);fn=fsam/2;
        W(1,1)=f1/fn; W(1,2)=f2/fn;
        [a b]=butter(2,W);
        ntime=0; %Initialize the next window
        nskip=ninici+1; %Initialize the next advancing
        %Move the analysis window through the seismic traces
        while ntime<num_time & nskip+tprime-1<=ndat
            ntime=ntime+1;
            u=dat(1,nskip:nskip+tprime-1)';
            %Filter the seismic traces with a band-pass
            u=filter(a,b,u);
            y1=u.*tukeywin(tprime);%tapered cosine windows
            %Apply the rms algorithm
            amp(ntime)=sqrt(sum((y1.*conj(y1))/tprime));
            %Increment temporal analysis variables
            nskip=nskip+round(time_step/delta_time);
            ssec=ssec+time_step;
            if ssec>=60
                ssec=ssec-60;
                imin=imin+1;
                if imin>=60
                    imin=imin-60;
                    ihr=ihr+1;
                    if ihr>=24
                        ihr=ihr-24;
                        iday=iday+1;
                    end
                end
            end
        end
        %Select the position of the maximum amplitude
        [~,imax]=max(amp);
    catch
        %Error assessment
        if f1>=round(fsam/2)|f2>=round(fsam/2);set(text1, 'String','Invalid analysis frequencies.');drawnow;return;end
    end
end

%==========================================================================
%%%%Plot of the FK-power values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting1(ii,ax)
    try
        nbins=str2num(get(h.nbins,'string'));timestep=hours(nbins);%time bin in hours
        cla(ax)%Clear the fk-power axes
        POW=cell2mat(fkpow);POW=POW(:,ii);%Select the ii-th block of the fk-power values
        Time=cell2mat(TIME_music);
        %Set the the limits and steps analysis
        time1=min(Time);time2=max(Time)+datenum(timestep);
        %Calculation and plot the 3D temporal histograms of  the fk-power values
        DistrHist3D(ax,Time,POW,time1,time2,timestep,0,1,0.1)
        ylim(ax,[0 1])
        xlim(ax,[datetime(time1,'ConvertFrom','datenum') datetime(max(Time),'ConvertFrom','datenum')])%Convert x-axis limits in datetime format
        ylabel(ax,'FK-POWER');
        legend(ax,['Signal-' num2str(slidval)])
    catch
        %Error assessment
        if size(POW,1)~=size(Time,1);set(text1, 'String','Error. Incoherent dimensions of the MUSIC array.');drawnow;return;end
    end
end

%==========================================================================
%%%%Plot of the BackAzimuth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting2(ii,ax)
    try
        nbins=str2num(get(h.nbins,'string'));timestep=hours(nbins);%time bin in hours
        cla(ax)%Clear the backazimuth axes
        BAZ=cell2mat(bazi);BAZ=BAZ(:,ii);%Select the ii-th block of the backazimuth
        Time=cell2mat(TIME_music);
        %Set the the limits and steps analysis
        time1=min(Time);time2=max(Time)+datenum(timestep);baz1=0;baz2=360;bazstep=10;
        %Calculation and plot the 3D temporal histograms of  the backazimuth
        DistrHist3D(ax,Time,BAZ,time1,time2,timestep,baz1,baz2,bazstep)
        ylim(ax,[0 360])
        xlim(ax,[datetime(time1,'ConvertFrom','datenum') datetime(max(Time),'ConvertFrom','datenum')])%Convert x-axis limits in datetime format
        ylabel(ax,'Backazimuth (degrees)');
        legend(ax,['Signal-' num2str(slidval)])
    catch
        %Error assessment
        if size(BAZ,1)~=size(Time,1);set(text1, 'String','Error. Incoherent dimensions of the MUSIC array.');drawnow;return;end
    end
end

%==========================================================================
%%%%Plot of the Ray parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting3(ii,ax)
    try
        nbins=str2num(get(h.nbins,'string'));timestep=hours(nbins);%time bin in hours
        cla(ax)%Clear the ray parameter axes
        SLOW=cell2mat(Slow);
        SLOW=SLOW(:,ii);%Clear the ray parameter axes
        Time=cell2mat(TIME_music);
        %Set the the limits and steps analysis
        time1=min(Time);time2=max(Time)+datenum(timestep);
        %Calculation and plot the 3D temporal histograms of the ray parameter
        DistrHist3D(ax,Time,SLOW,time1,time2,timestep,0,10,0.1)
        ylim(ax,[0 5])
        xlim(ax,[datetime(time1,'ConvertFrom','datenum') datetime(max(Time),'ConvertFrom','datenum')])%Convert x-axis limits in datetime format
        ylabel(ax,'RayParameter (s/km)');
        legend(ax,['Signal-' num2str(slidval)])
    catch
        %Error assessment
        if size(SLOW,1)~=size(Time,1);set(text1, 'String','Error. Incoherent dimensions of the MUSIC array.');drawnow;return;end
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
%%%%Refresh the MUSIC axes through the pre button%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loca1(~,~)
   %Error control
   if isempty(cell2mat(bazi));set(text1, 'String','Error. MUSIC matrix is empty.');return;end
   %Decrease the slidval value
   slidval=slidval-1;if slidval<=0;slidval=1;end
   %Plot the new block of the fk-power axes
   ploting1(slidval,axes1);
   %Plot the new block of the backazimuth
   ploting2(slidval,axes2);
   %Plot the new block of the ray parameter
   ploting3(slidval,axes3);
        
end
%==========================================================================
%%%%Refresh the MUSIC axes through the next button%%%%%%%%%%%%%%%%%%%%%%%%%
function loca2(~,~)
   %Error control
   if isempty(cell2mat(bazi));set(text1, 'String','Error. MUSIC matrix is empty.');return;end
   %Increase the slidval value
   slidval=slidval+1;if slidval>nsig;slidval=nsig;end
   %Plot the new block of the fk-power axes
   ploting1(slidval,axes1);
   %Plot the new block of the backazimuth
   ploting2(slidval,axes2);
   %Plot the new block of the ray parameter
   ploting3(slidval,axes3);        
end

%==========================================================================
%%%%Set the limits of the y-axis on the fk-power axes%%%%%%%%%%%%%%%%%%%%%%
function setlimits1(~,~)
    try
        %Manual settings
        %Acquisition limits
        pow1_lim=str2num(get(h.pow1_lim,'string'));
        pow2_lim=str2num(get(h.pow2_lim,'string'));  
        ylim(axes1,[pow1_lim pow2_lim])
    catch
        %Default settings
        ylim(axes1,[0 1])
        %Reset fk power cross correlation limits
        set(h.pow1_lim,'string','0');
        set(h.pow2_lim,'string','1');
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
%%%%Set the limits of the x-axis on the MUSIC axes%%%%%%%%%%%%%%%%%%%%%%%%%
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
        if ~isempty(cell2mat(TIME_music))
            Time=datetime(cell2mat(TIME_music),'ConvertFrom','datenum');%Convert in datetime format
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
    freq0=str2num(get(h.freq0,'string'));
    nfreq=str2num(get(h.nfreq,'string'));
    m=str2num(get(h.m,'string'));
    lwin=str2num(get(h.lwin,'string'));
    advance=str2num(get(h.advance,'string'));
    NSig1=str2num(get(h.nsig,'string'));
    eigperc=str2num(get(h.eigperc,'string'));
    R=str2num(get(h.R,'string'));
    nbins=str2num(get(h.nbins,'string'));
    vel=str2num(get(h.vel,'string'));
    
    %Control and eventually reset parameters
    if (size(freq0,1)==0);set(h.freq0,'string',params.freq0);end
    if (size(nfreq,1)==0);set(h.nfreq,'string',params.nfreq);end
    if (size(m,1)==0);set(h.m,'string',params.m);end
    if (size(lwin,1)==0);set(h.lwin,'string',params.lwin);end
    if (size(advance,1)==0);set(h.advance,'string',params.advance);end
    if (size(m,1)==0);set(h.m,'string',params.m);end
    if (size(NSig1,1)==0);set(h.nsig,'string',params.nsig);end
    if (size(eigperc,1)==0);set(h.eigperc,'string','95');end
    if (size(R,1)==0);set(h.R,'string',params.R);end
    if (size(vel,1)==0);set(h.vel,'string','1.6');end
    if (size(nbins,1)==0);set(h.nbins,'string','100');end
end


%==========================================================================
%%%%Save MUSIC plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save1(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    try
        nbins=str2num(get(h.nbins,'string'));timestep=hours(nbins);%time bin in hours
        ii=slidval;
        POW=cell2mat(fkpow);
        BAZ=cell2mat(bazi);
        SLOW=cell2mat(Slow);
        Time=cell2mat(TIME_music);
        POW=POW(:,ii);%Select the ii-th block of the fk-power
        BAZ=BAZ(:,ii);%Select the ii-th block of the backazimuth
        SLOW=SLOW(:,ii);%Select the ii-th block of the ray parameter
        %Set the the limits and steps analysis
        time1=min(Time);time2=max(Time)+datenum(seconds(timestep));
        new_figure=figure('Name','MUSIC plot','NumberTitle','off')
        %Plot the values of the fk-power
        subplot(3,1,1)
        %Calculation and plot the 3D temporal histograms of the fk-power values
        NDistrHist3D(Time,POW,time1,time2,timestep,0,1,0.1)
        ylim([0 1])
        xlim([datetime(time1,'ConvertFrom','datenum') datetime(max(Time),'ConvertFrom','datenum')])%Convert x-axis limits in datetime format
        ylabel('FK-POWER');
        legend(['Signal-' num2str(slidval)])
        %Plot the values of the backazimuth
        subplot(3,1,2)
        %Calculation and plot the 3D temporal histograms of the backazimuth
        NDistrHist3D(Time,BAZ,time1,time2,timestep,0,360,10)
        ylim([0 360])
        xlim([datetime(time1,'ConvertFrom','datenum') datetime(max(Time),'ConvertFrom','datenum')])%Convert x-axis limits in datetime format
        ylabel('BackAzimuth (degrees)');
        legend(['Signal-' num2str(slidval)])
        %Plot the values of the rayparameter
        subplot(3,1,3)
        %Calculation and plot the 3D temporal histograms of the rayparameter
        NDistrHist3D(Time,SLOW,time1,time2,timestep,0,10,0.1)
        ylim([0 5])
        xlim([datetime(time1,'ConvertFrom','datenum') datetime(max(Time),'ConvertFrom','datenum')])%Convert x-axis limits in datetime format
        ylabel('RayParameter (s/km)');
        legend(['Signal-' num2str(slidval)])
    catch
        %Error assessment
        if isempty(cell2mat(Time)) | isempty (cell2mat(BAZ));set(text1, 'String','Error. MUSIC matrix is empty.');drawnow;return;end
        if size(BAZ,1)~=size(Time,1);set(text1, 'String','Error. Incoherent dimensions of the MUSIC array.');drawnow;return;end
    end
end

%==========================================================================
%%%%Plot the 3D Multivariate analysis Ray parameter vs. Backazimuth%%%%%%%%
function save2(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    try
        ii=slidval;
        BAZ=cell2mat(bazi);
        Ray=cell2mat(Slow);
        BAZ=BAZ(:,ii);%Select the ii-th block of the backazimuth
        Ray=Ray(:,ii);%Select the ii-th block of the ray parameter
        %Set the the limits and steps analysis
        az1=0;az2=360;azstep=10;ray1=0;ray2=10;raystep=0.1;
        new_figure=figure('Name','Ray Parameter vs. BackAzimuth','NumberTitle','off');
        %Calculation and plot the 3D Multivariate analysis
        NewHist3D(Ray,BAZ,ray1,ray2,raystep,az1,az2,azstep)
        set(gca,'ytick', [0:36:360], 'ylim', [0 360])
        set(gca,'xtick', [0:0.5:5], 'xlim', [0 5])
        ylabel('BackAzimuth (degress)');
        xlabel('RayParameter (s/km)');
        legend(['Signal-' num2str(slidval)])  
    catch
        %Error assessment
        if isempty(cell2mat(Ray)) | isempty (cell2mat(BAZ));set(text1, 'String','Error. MUSIC matrix is empty.');drawnow;return;end
        if size(BAZ,1)~=size(Ray,1);set(text1, 'String','Error. Incoherent dimensions of the MUSIC array.');drawnow;return;end
    end
end

%==========================================================================
%%%%Plot the 3D Multivariate analysis of the horizontal slowness%%%%%%%%%%%
function save3(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    try
        ii=slidval;
        sx=cell2mat(Sxmax);
        sy=cell2mat(Symax);
        sx=sx(:,ii);sy=sy(:,ii);%Select the ii-th block of the horizontal slowness
        new_figure=figure('Name','Slowness Grid','NumberTitle','off');
        %Calculation and plot the 3D Multivariate analysis
        NewHist3D(sx,sy,-params.pmax,params.pmax,0.1,-params.pmax,params.pmax,0.1)
        set(gca,'ytick', [-params.pmax:0.8:params.pmax], 'ylim', [-params.pmax params.pmax])
        set(gca,'xtick', [-params.pmax:0.8:params.pmax], 'xlim', [-params.pmax params.pmax])
        ylabel('Sy (s/km)');
        xlabel('Sx (s/km)');
        legend(['Signal-' num2str(slidval)]) 
    catch
        %Error assessment
        if isempty(cell2mat(sx)) | isempty (cell2mat(sy));set(text1, 'String','Error. MUSIC matrix is empty.');drawnow;return;end
        if size(sx,1)~=size(sy,1);set(text1, 'String','Error. Incoherent dimensions of the MUSIC array.');drawnow;return;end     
    end
end

%==========================================================================
%%%%Plot the histograms of the Incidence and Backazimuth%%%%%%%%%%%%%%%%%%%
function save4(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    try
        ii=slidval;
        INC=cell2mat(Inc);
        BAZ=cell2mat(bazi);
        INC=INC(:,ii);%Select the ii-th block of the incidence
        BAZ=BAZ(:,ii);%Select the ii-th block of the backazimuth
        %Incidence bins
        bins=0:10:90;
        %Backazimuth bins
        fin=(2*pi*BAZ)./360;
        bins2=[0:5*((2*pi)/360):355*((2*pi)/360)];
        new_figure=figure('Name','Incidence vs. Back Azimuth','NumberTitle','off');
        %Plot the histogram of the incidence
        subplot(2,1,1)
        newh=histogram(INC,bins,'FaceColor','b');
        newh.Normalization = 'probability';
        xlim([0 90])
        ylabel('Probability')
        legend(['Signal-' num2str(slidval)])
        grid on; grid minor;
        %Plot the polar histogram of the backazimuth
        subplot(2,1,2)
        newh2=polarhistogram(fin,bins2);
        newh2.Normalization = 'probability';
        title(['Backazimuth (degrees)'])
        set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','clockwise');
        legend(['Signal-' num2str(slidval)])  
    catch
        %Error assessment
        if isempty(cell2mat(INC)) | isempty (cell2mat(BAZ));set(text1, 'String','Error. MUSIC matrix is empty.');drawnow;return;end
        if size(BAZ,1)~=size(INC,1);set(text1, 'String','Error. Incoherent dimensions of the MUSIC array.');drawnow;return;end  
    end
end

%==========================================================================
%%%%Save the results of the MUSIC analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_music(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Validation of the results
    if isempty(params.save);set(text1, 'String', sprintf('Error. No output directory.'));drawnow;return;end  
    if isempty(cell2mat(bazi));set(text1, 'String', sprintf('Error. MUSIC matrix is empty.'));drawnow;return;end
    if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
    %Loop through the sources number
    for jj=1:size(bazi,1)
        %Create the Output filename
        name=ff{jj};[~,name,~] = fileparts(name);
        switch params.comp_chan
            case 'Chan'
                name=name(end-14:end);
            case 'Comp'  
                name=name(end-18:end);
        end
        Baz=bazi{jj};Ray=Slow{jj};Pow=fkpow{jj};ssx=Sxmax{jj};ssy=Symax{jj};incc=Inc{jj};Time=TIME_music{jj};
        %Loop through the number of Input files (ff)
        for kk=1:size(Baz,2)
            time=Time;baz=Baz(:,kk);rayp=Ray(:,kk);pow=Pow(:,kk);sx=ssx(:,kk);sy=ssy(:,kk);inc=incc(:,kk);
            %Compute the juliand day
            jday=day(datetime(time(1),'ConvertFrom','datenum'),'dayofyear');
            if jday<10;jday=['00' num2str(jday)];elseif jday>=10 & jday<100;jday=['0' num2str(jday)];elseif jday>=100;jday=[num2str(jday)];end
            %Create the Output directory
            Ffolder=['Sign-' num2str(kk)];
            pathout=[params.save '/MUSIC/' Ffolder '/' jday '/'];if ~exist(pathout, 'dir');mkdir(pathout);end 
            output=[pathout name params.outputformat];
            %Create a table of the results
            data=table(time,baz,rayp,inc,sx,sy,pow);data.Properties.VariableNames={'Time','BackAzimuth', 'RayParameter','Incidence','Sx','Sy','FkPower'};
            data.Properties.VariableUnits = {'datenum','degrees','s/km ','degrees','s/km ','s/km ',' '};
            %Save the results
            switch params.outputformat
                case '.mat'
                    save (output,'data');
                case '.txt'
                    writetable(data,output)
            end
        end
    end
    %Command message
    set(text1, 'String', sprintf(['MUSIC' params.outputformat '-Observed data are saved.']));
    drawnow
end
end 