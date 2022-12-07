%%%%%%%%%%%%% Create the main Data structures%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params: user parameters
% avfig: Create Data Structure figure
% panel_color: colors of the main panel 
% fig_color: colors of the figure
% pan1: control right side panel 
% text1: Window message
% h,h2: UI control parameters
%%%%%%%%%%%h.variables/h2.variables/params.variables/variables%%%%%%%%%%%%%
% rad,Rad,RRad: radio button within a button group 
% rad1, rad2, rad3, rad4, Rad1, Rad2, RRad1, RRad2, RRad3: button groups 
%-------------------------------- Prexisting files settings
% data, sacin, xlm: the Input folder of the prexisting files 
% save,savein, savename: the Output folder of the prexisting files 
% format_file, format: the format (.sac/.mseed/.cube) of the Input files
% sta_ref, staz_reff,sta_coord: the name of the reference station 
% component, comp: the name of the component/channel
% cchan,comp_chan: the station system used (component/channel)
% lat_coord, lon_coord, ele_coord: coordinates in wgs84 system
% pols, poles: poles of the transfer function
% zers, zeros: zeros of the transfer function
% k: calibration coefficient (rad/s)
% vvel: coefficient for conversion from Volts to m/s (V s/m)
% bmw: coefficient for conversion from counts to Volts (V/counts)
%-------------------------------- Iris files settings
% irisout, data2: the Output folder of the Iris files
% sta: the name of the reference station 
% net: network of reference
% loc: localization of reference
% channe: channel of reference
% t1,t2: temporal limits of reference
%-------------------------------- General settings
% search_sac, date_search: term of research file 
% pathname/pathnam/pathaname2, fname: Input folder and list of filenames
% ff: list of files with their path
% list2: list of seismic traces directories
% mytrace: seismic trace structure
% dt: sampling interval
% fs: sampling rate
% nyq: nyquist number
% z: amplitude vector of the seismic trace
% nt: number of samples
% startTime, endTime: temporal limits of the traces 
% jday, hh, mi,ss: temporal variables 

function CreateDataStruct()
%% Declaration global parameters
global h2 dispval fval pathname pathnam pathname2 data ff comp text1 data2 savename2 params xlm pathsta ValMode comps_ref
ff= cell(1);
dispval=0;fval=0;fs=[];pathname=[];pathnam=[];pathname2=[];pathsta=[];
dt=[];icons=[]; mytrace=[]; comp_chan=[]; format=[];ext=[];component=[];cchan=[];format_file=[];
load('icons.mat')
%Set Input file directories
data='';
data2='';
%Set Output file directories
savename='';
xlm='';
%Set the station system used
comp=params.component;%%Set Component
comp_chan=params.comp_chan;
format='Sac ';
comps_ref={'Z','N','E','F'};
ValMode=2;
%% Create Data Structures figure
%%Figure
fig_color= [.94 .94 .94];
panel_color= [.97 .97 .97];

avfig= figure('Name', 'CreateDataStructures', 'position',get(0,'screensize'), ...
    'NumberTitle','off','units','normalized', 'color', fig_color,'MenuBar','none');%%Create the main figure

%% Control right side panel            
pan1= uipanel(avfig,'visible','on','Position',[.01 .01 .98 .98],...
    'BackgroundColor', panel_color);

%Window message to user
text1= uicontrol('Style','edit','FontName','Arial','Fontsize',10,'Fontweight','bold',...
                'String','Window Message','Max',5,'Min',0,'HorizontalAlignment','left',...
                'BackgroundColor', 'w','Enable','inactive',...
                'units','normalized','Position',[.165 .05 .72 .05]);
  

%% Text boxes
%-------------------------------- Prexisting files setting
uicontrol(pan1,'Style','text', 'String','Prexisting file',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.05 .92 .8 .07],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Input',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.08 .9 .2 .05],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Output',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.08 .855 .2 .05],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Station stats',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.08 .805 .2 .05],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','XML file',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.08 .755 .2 .05],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Staz',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.08 .74 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Lat (°)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.25 .74 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Lon (°)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.37 .74 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Ele(m)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.49 .74 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','k (rad/s)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.60 .74 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','C2V (V/counts)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.715 .74 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','S (V s/m)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.86 .74 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Poles (rad/s)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.62 .655 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Zeros (rad/s)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.8 .655 .2 .02],...
    'BackgroundColor',panel_color);

%-------------------------------- Iris files setting
uicontrol(pan1,'Style','text', 'String','IRIS file',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.05 .48 .8 .07],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Output',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.08 .44 .2 .05],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Net',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.08 .38 .2 .05],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Staz',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.08 .36 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Loc',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.08 .30 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Chan',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.08 .24 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','t1',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.08 .195 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','t2',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.08 .145 .2 .02],...
    'BackgroundColor',panel_color);

%% Editable text 
%-------------------------------- Prexisting files setting
h2.sacin= uicontrol(pan1,'Style','edit', 'String',data,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.16 .925 .73 .04],...
    'BackgroundColor','w');
h2.sacinicon= uicontrol(pan1,'Style','pushbutton', 'String','',...
    'Cdata',icons.folder,'Tag','sacin',...
    'Units','normalized','Position',[.9 .925 .06 .04],...
    'callback',@folder_cbk);
h2.savein= uicontrol(pan1,'Style','edit', 'String',savename,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.16 .875 .73 .04],...
    'BackgroundColor','w');
h2.saveinicon= uicontrol(pan1,'Style','pushbutton', 'String','',...
    'Cdata',icons.folder,'Tag','savein',...
    'Units','normalized','Position',[.9 .875 .06 .04],...
    'callback',@folder_cbk);
h2.pathsta= uicontrol(pan1,'Style','edit', 'String',pathsta,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.16 .825 .73 .04],...
    'BackgroundColor','w');
h2.pathstaicon= uicontrol(pan1,'Style','pushbutton', 'String','',...
    'Cdata',icons.folder,'Tag','pathsta',...
    'Units','normalized','Position',[.9 .825 .06 .04],...
    'callback',@folder_cbk);
h2.xlmicon= uicontrol(pan1,'Style','pushbutton', 'String','',...
    'Cdata',icons.file,'Tag','xlm',...
    'Units','normalized','Position',[.9 .775 .06 .04],...
    'Enable','off','callback',@folder_cbk);
h2.xlm= uicontrol(pan1,'Style','edit', 'String',xlm,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.16 .775 .73 .04],...
    'BackgroundColor','w','Enable','off');
h2.staz_reff= uicontrol(pan1,'Style','edit', 'String',params.sta_ref,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.16 .73 .06 .04],...
    'BackgroundColor','w','callback',@setlimits1);
h2.lat_coord= uicontrol(pan1,'Style','edit', 'String',[],...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.28 .73 .06 .04],...
    'BackgroundColor','w','Enable','off','callback',@setlimits2);
h2.lon_coord= uicontrol(pan1,'Style','edit', 'String',[],...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.4 .73 .06 .04],...
    'BackgroundColor','w','Enable','off','callback',@setlimits2);
h2.ele_coord= uicontrol(pan1,'Style','edit', 'String',[],...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.52 .73 .06 .04],...
    'BackgroundColor','w','Enable','off','callback',@setlimits2);
h2.k= uicontrol(pan1,'Style','edit', 'String',[],...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.64 .73 .06 .04],...
    'BackgroundColor','w','Enable','off','callback',@setlimits2);
h2.bmw= uicontrol(pan1,'Style','edit', 'String',[],...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.78 .73 .06 .04],...
    'BackgroundColor','w','Enable','off','callback',@setlimits2);
h2.vvel= uicontrol(pan1,'Style','edit', 'String',[],...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.90 .73 .06 .04],...
    'BackgroundColor','w','Enable','off','callback',@setlimits2);
h2.pols= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left','Max',2,...
    'Units','normalized','Position',[.68 .605 .1 .1],...
    'BackgroundColor','w','Enable','off');
h2.zers= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left','Max',2,...
    'Units','normalized','Position',[.86 .605 .1 .1],...
    'BackgroundColor','w','Enable','off');

%-------------------------------- Iris files setting
h2.irisout= uicontrol(pan1,'Style','edit', 'String',data2,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.16 .455 .6 .04],...
    'BackgroundColor','w');
h2.irisicon= uicontrol(pan1,'Style','pushbutton', 'String','',...
    'Cdata',icons.folder,'Tag','irisout',...
    'Units','normalized','Position',[.77 .455 .06 .04],...
    'callback',@folder_cbk);
h2.net= uicontrol(pan1,'Style','edit', 'String','IU',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.16 .4 .3 .04],...
    'BackgroundColor','w','callback',@setlimits2);
h2.sta= uicontrol(pan1,'Style','edit', 'String','ANMO',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.16 .345 .3 .04],...
    'BackgroundColor','w','callback',@setlimits2);
h2.loc= uicontrol(pan1,'Style','edit', 'String','00',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.16 .290 .3 .04],...
    'BackgroundColor','w','callback',@setlimits2);
h2.channe= uicontrol(pan1,'Style','edit', 'String','BHZ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.16 .235 .3 .04],...
    'BackgroundColor','w','callback',@setlimits2);
h2.t1= uicontrol(pan1,'Style','edit', 'String','2010-02-27 06:00:00',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.16 .180 .3 .04],...
    'BackgroundColor','w','callback',@setlimits2);
h2.t2= uicontrol(pan1,'Style','edit', 'String','2010-02-27 07:00:00',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.16 .125 .3 .04],...
    'BackgroundColor','w','callback',@setlimits2);
h2.dip=[];
h2.az=[];
h2.depth=[];
h2.units=[];

%% Add Radio buttons
%-------------------------------- Prexisting files setting
% Editable text
uicontrol(pan1,'Style','text', 'String','Comp/Chan',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.08 .67 .4 .05],...
    'BackgroundColor',panel_color);

% Station system selection radio buttons
h.Rad= uibuttongroup(pan1,'units','normalized','BackgroundColor',panel_color,...
    'bordertype','none','Position',[.13 .68 .24 .05]);
set(h.Rad,'SelectionChangeFcn',@radcbk3);
h.Rad1 = uicontrol( h.Rad, 'Style','Radio','String','Comp',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.15 0 .4 1],'HandleVisibility','off');
h.Rad2 = uicontrol( h.Rad, 'Style','Radio','String','Chan',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.3 0 .4 1],'HandleVisibility','off');

% Set the station system of reference
if strcmp(comp_chan,'Comp')
    set (h.Rad1,'value',1);
    set (h.Rad2,'value',0);

elseif strcmp(comp_chan,'Chan')
    set (h.Rad1,'value',0);
    set (h.Rad2,'value',1);
end

% Component selection radio buttons
h.rad= uibuttongroup(pan1,'units','normalized','BackgroundColor',panel_color,...
    'bordertype','none','Position',[.27 .68 .24 .05]);
set(h.rad,'SelectionChangeFcn',@radcbk);
h.rad1 = uicontrol( h.rad, 'Style','Radio','String','Z',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.0 0 .4 1],'HandleVisibility','off');
h.rad2 = uicontrol( h.rad, 'Style','Radio','String','N',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.1 0 .4 1],'HandleVisibility','off');
h.rad3 = uicontrol( h.rad, 'Style','Radio','String','E',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.2 0 .4 1],'HandleVisibility','off');
h.rad4 = uicontrol( h.rad, 'Style','Radio','String','F',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.3 0 .4 1],'HandleVisibility','off');

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

% Format file selection radio buttons
h.RRad= uibuttongroup(pan1,'units','normalized','BackgroundColor',panel_color,...
    'bordertype','none','Position',[.41 .68 .24 .05]);
set(h.RRad,'SelectionChangeFcn',@radcbk2);
h.RRad1 = uicontrol( h.RRad, 'Style','Radio','String','Sac ',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.0 0 .4 1],'HandleVisibility','off');
h.RRad2 = uicontrol( h.RRad, 'Style','Radio','String','Seed',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.15 0 .4 1],'HandleVisibility','off');
h.RRad3 = uicontrol( h.RRad, 'Style','Radio','String','Cube',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.3 0 .4 1],'HandleVisibility','off');

%Set the format of reference 
if strcmp(format,'Sac ')
    set (h.RRad1,'value',1);
    set (h.RRad2,'value',0);
    set (h.RRad3,'value',0);

elseif strcmp(format,'Seed')
    set (h.RRad1,'value',0);
    set (h.RRad2,'value',1);
    set (h.RRad3,'value',0);

elseif strcmp(format,'Cube')
    set (h.RRad1,'value',0);
    set (h.RRad2,'value',0);
    set (h.RRad3,'value',1);
end
%% Buttons
%-------------------------------- Prexisting files setting
% Open files button
uicontrol('style','pushbutton', 'string','Open files','units','normalized',...
    'position',[.16 .56 .1 .06], 'callback',@openfnc);
% Open folder button
uicontrol('style','pushbutton', 'string','Open folder','units','normalized',...
    'position',[.285 .56 .1 .06], 'callback',@openfolder);
% Save button
uicontrol('style','pushbutton', 'string','Save','units','normalized',...
    'position',[.41 .56 .1 .06], 'callback',@save1); 

% Single/Multiple data preparation checkbox
h2.delim= uicontrol(pan1,'Style','popup', ...
    'String',{'Manual','Automatic'},'Value',2,...
    'HorizontalAlignment','left','Enable','on',...
    'Units','normalized','Position',[.16 .895 .3 .1],...
    'BackgroundColor','w','callback',@setMode);

%-------------------------------- Iris files setting
uicontrol('style','pushbutton', 'string','Save','units','normalized',...
    'position',[.55 .26 .1 .06], 'callback',@save2);

%% Other functions (nested)
%%%%Return the list of files manually selected%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function openfnc(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    staz_ref=deblank(get(h2.staz_reff,'string'));%Get the station name of reference
    data=get(h2.sacin,'string');%Get the Input directory
    tt= 0;pathnam=[];pathname2=[];
    ff=cell(1);%Initialize the list of files
    i=0;
    component=comp;
    cchan=comp_chan;
    %Define the search term on the basis of the format file
    switch format
        case 'Sac '
            search_sac=strcat('*',staz_ref,'*',component,'*.sac');
        case 'Seed'
            search_sac=strcat('*',staz_ref,'*',component,'*');
        case 'Cube'
            search_sac=strcat('*',staz_ref,'*');         
    end
    format_file=format;
    [fname, pathname] = uigetfile({search_sac;'*.*'},'File Selector',...
        data,'MultiSelect','on');%Selection of the files      
    
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
    text1.String=[];set(text1);drawnow;%Clear the Window message
    pathname=[];
    data=get(h2.sacin,'string');%Get the Input directory
    pathnam= uigetdir(data);
    pathnam=[pathnam '/'];
    staz_ref=deblank(get(h2.staz_reff,'string'));%Get the station name of reference 
    component=comp;
    cchan=comp_chan;
    i=0;
    %Define the search term on the basis of the format file
    switch format
        case 'Sac '
            search_sac=strcat('*',staz_ref,'*', component,'*.sac');
        case 'Seed'
            search_sac=strcat('*',staz_ref,'*', component,'*');
        case 'Cube'
            search_sac=strcat('*',staz_ref,'*');         
    end
    format_file=format;
    if pathnam~=0
    % Get all files in the directory and subdirectories
    [ff,pathname2]= getAllFile(pathnam,search_sac);
    ff=ff';pathname2=pathname2';
    
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
    if ~isempty(h2.lat_coord.String);set(h2.lat_coord,'string',[]);end
     if ~isempty(h2.lon_coord.String);set(h2.lon_coord,'string',[]);end
      if ~isempty(h2.ele_coord.String);set(h2.ele_coord,'string',[]);end
     if ~isempty(h2.k.String);set(h2.k,'string',[]);end
     if ~isempty(h2.bmw.String);set(h2.bmw,'string',[]);end
      if ~isempty(h2.vvel.String);set(h2.vvel,'string',[]);end
      if ~isempty(h2.xlm.String);set(h2.xlm,'string',[]);end
      if ~isempty(h2.pols.String);set(h2.pols,'string',[]);end
      if ~isempty(h2.zers.String);set(h2.zers,'string',[]);end
      ff=cell(1);
    
end

%==========================================================================
%%%%%%%%%%%%%%%%Return the format of the Input files%%%%%%%%%%%%%%%%%%%%%%%
function radcbk2(~,eventdata)
    
    format = get(eventdata.NewValue,'String');   
    
end

%==========================================================================
%%%%Return the type of the station system%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radcbk3(~,eventdata)
    
   comp_chan = get(eventdata.NewValue,'String');   
    
end

%==========================================================================
%%%%Control and reset automatically station name%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setlimits1(~,~)
    %Acquisition station name
    staz_ref=deblank(get(h2.staz_reff,'string'));
    ff=cell(1);%Initialize the list of files
    %Control and eventually reset parameter
    if ~isempty(h2.lat_coord.String);set(h2.lat_coord,'string',[]);end
     if ~isempty(h2.lon_coord.String);set(h2.lon_coord,'string',[]);end
      if ~isempty(h2.ele_coord.String);set(h2.ele_coord,'string',[]);end
     if ~isempty(h2.k.String);set(h2.k,'string',[]);end
     if ~isempty(h2.bmw.String);set(h2.bmw,'string',[]);end
      if ~isempty(h2.vvel.String);set(h2.vvel,'string',[]);end
      if ~isempty(h2.xlm.String);set(h2.xlm,'string',[]);end
      if ~isempty(h2.pols.String);set(h2.pols,'string',[]);end
      if ~isempty(h2.zers.String);set(h2.zers,'string',[]);end
      ff=cell(1);
    if (size(staz_ref,2)==0);set(h2.staz_reff,'string',params.sta_ref);end
end

%==========================================================================
%%%%Control and reset automatically parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setlimits2(~,~)
    %Acquisition parameters
    lat_coord=get(h2.lat_coord,'string');
    lon_coord=get(h2.lon_coord,'string');
    ele_coord=get(h2.ele_coord,'string');
    sta=deblank(get(h2.sta,'string'));
    net=deblank(get(h2.net,'string'));
    loc=deblank(get(h2.loc,'string'));
    channe=deblank(get(h2.channe,'string'));
    t1=deblank(get(h2.t1,'string'));
    t2=deblank(get(h2.t2,'string'));
    
    %Control and eventually reset parameters
    if (size(lat_coord,2)==0);set(h2. lat_coord,'string',[]);end
    if (size(lon_coord,2)==0);set(h2. lon_coord,'string',[]);end
    if (size(ele_coord,2)==0);set(h2. ele_coord,'string',[]);end
    if (size(sta,2)==0);set(h2.sta,'string','ANMO');end
    if (size(net,2)==0);set(h2.net,'string','IU');end
    if (size(loc,2)==0);set(h2.loc,'string','00');end
    if (size(channe,2)==0);set(h2.channe,'string','BHZ');end
    if (size(t1,2)==0);set(h2.t1,'string','2010-02-27 06:00:00');end
    if (size(t2,2)==0);set(h2.t2,'string','2010-02-27 07:00:00');end
end


%==========================================================================
%%%%%%%Set Manual/Automatic preparation of the datatset %%%%%%%%%%%%%%%%%%%
function setMode(h, ~)
    if get(h,'Value')==1
        ValMode=1;
        %Manual Mode
        set(h2.pathsta,'Enable','off');
         set(h2.pathstaicon,'Enable','off');
         set(h2.lat_coord,'Enable','on');
         set(h2.lon_coord,'Enable','on');
         set(h2.ele_coord,'Enable','on');
         set(h2.k,'Enable','on');
         set(h2.bmw,'Enable','on');
         set(h2.vvel,'Enable','on');
         set(h2.pols,'Enable','on');
         set(h2.zers,'Enable','on');
         set(h2.xlm,'Enable','on');
         set(h2.xlmicon,'Enable','on');
         if ~isempty(h2.pathsta.String);set(h2.pathsta,'string',[]);end
         ff=cell(1);
    elseif get(h,'Value')==2
        ValMode=2;
        %Automatic Mode
        set(h2.pathsta,'Enable','on'); 
        set(h2.pathstaicon,'Enable','on');
        set(h2.lat_coord,'Enable','off');
        set(h2.lon_coord,'Enable','off');
        set(h2.ele_coord,'Enable','off');
        set(h2.k,'Enable','off');
        set(h2.bmw,'Enable','off');
        set(h2.vvel,'Enable','off');
        set(h2.pols,'Enable','off');
        set(h2.zers,'Enable','off');
        set(h2.xlm,'Enable','off');
        set(h2.xlmicon,'Enable','off');
        if ~isempty(h2.lat_coord.String);set(h2.lat_coord,'string',[]);end
        if ~isempty(h2.lon_coord.String);set(h2.lon_coord,'string',[]);end
        if ~isempty(h2.ele_coord.String);set(h2.ele_coord,'string',[]);end
        if ~isempty(h2.k.String);set(h2.k,'string',[]);end
        if ~isempty(h2.bmw.String);set(h2.bmw,'string',[]);end
        if ~isempty(h2.vvel.String);set(h2.vvel,'string',[]);end
        if ~isempty(h2.xlm.String);set(h2.xlm,'string',[]);end
        if ~isempty(h2.pols.String);set(h2.pols,'string',[]);end
         if ~isempty(h2.zers.String);set(h2.zers,'string',[]);end
        ff=cell(1);
    end  
end

%==========================================================================
%%%%Convert and Save the prexisting files in Matlab format%%%%%%%%%%%%%%%%%
function save1(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    try
        %Error assessment
        if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end 
        if ValMode==1
        %%%%%%%%Manual Mode%%%%%%%%%%%%%%%%%%%%%%%
        savename=get(h2.savein,'string');%Get the Output directory
        %Error control
        if isempty(savename);set(text1, 'String', sprintf('Error. No output directory.'));drawnow;return;end
        %Command message
        textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
        set(text1, 'String', textLabel); 
        drawnow
        nsize=length(ff);%Number of files selected
        %Acquisition and validation of the parameters
        staz_ref=deblank(get(h2.staz_reff,'string'));
        lat_coord=str2num(get(h2.lat_coord,'string'));
        lon_coord=str2num(get(h2.lon_coord,'string'));
        ele_coord=str2num(get(h2.ele_coord,'string'));
        dip=h2.dip;
        az=h2.az;
        depth=h2.depth;
        units=h2.units;
        k=str2num(get(h2.k,'string'));if k<=0;set(text1, 'String','Error. Invalid normalization factor value.');drawnow;return;end
        vvel=str2num(get(h2.vvel,'string'));if vvel<=0;set(text1, 'String','Error. Invalid sensivity value.');drawnow;return;end
        bmw=str2num(get(h2.bmw,'string'));if bmw<=0;set(text1, 'String','Error. Invalid  conversion counts to volt value.');drawnow;return;end
        zeros=str2num(get(h2.zers,'string'));
        poles=str2num(get(h2.pols,'string'));
        %Redefine the station name of reference
        switch cchan
            case 'Comp'
                sta_coord=staz_ref;
            case 'Chan'
                sta_coord=[staz_ref '.BH' component];
        end
        %Load the initial signal structure
        load (strcat('mytrace.mat'))
        %Execute one of three groups of routines (.sac,.mseed,.cube files)
        switch format_file
            case 'Sac '
                %Loop through the number of files selected
                for kk=1:nsize
                    %Read and load the kk-th file .sac 
                    filename=ff{kk};
                    T=rsac(filename); 
                    %Return the amplitude, time and info vectors
                    z=T(:,2);dt=T(1,3); fs=(1/dt);yy=T(71,3);jday=T(72,3);hh=T(73,3);mi=T(74,3);ss=T(75,3)+T(76,3)/1000;nt=length(z);
                     %Conversion from counts to Volts
                        if ~isnan(bmw) & ~isempty(bmw) 
                            z = z.*bmw;
                        end
                        %Conversion from Volts to m/s
                        if ~isnan(vvel) & ~isempty(vvel) 
                            z = z./vvel;
                        end    
                    %Define temporal variables and names
                        timeref=datenum(yy,0,jday,hh,mi,ss);
                        timeref=datetime(timeref,'ConvertFrom','datenum');
%                        timeref=datetime(timeref,'ConvertFrom','datenum');
                        timeref.Format='dd-MMM-yyyy HH:mm:SS';
                        jday=day(timeref,'dayofyear');
                        timecontrol1=datetime(datenum(yy,0,jday,hh,0,0),'ConvertFrom','datenum');
                        timecontrol1.Format='dd-MMM-yyyy HH:mm:SS';
                        timecontrol2=datetime(datenum(yy,0,jday,hh+1,0,0),'ConvertFrom','datenum');
                        timecontrol2.Format='dd-MMM-yyyy HH:mm:SS';
                        timerefend=timeref+seconds((length(z)/round(fs)));
                        timerefend.Format='dd-MMM-yyyy HH:mm:SS';
                        TimeThr=timecontrol1+seconds(1/round(fs)):seconds(1/round(fs)):timecontrol2;
                        idxtime1=find(TimeThr==timeref);
                        idxtime2=find(TimeThr==timerefend);
                        %Fill temporal gaps to create hourly files 
                        if nt~=3600*round(fs)
                            if timeref>timecontrol1
                                gapp=idxtime1;
                                z=[ones(gapp,1).*0;z];
                            end
                            if timerefend<timecontrol2
                                gapp=3600*round(fs)-idxtime2;
                                 z=[z;ones(gapp,1).*0];
                            end
                            timeref=timecontrol1;
                            mi=0;
                            ss=0;
                        end
                        startTime=datenum(timeref);
                        endTime=datenum(timeref+seconds((length(z)/round(fs)))); 
%                         endTime=t(1)+(length(d)/round(fs))/(24*3600)+(1/fs)/(24*3600);
                        if hh<10;hh=['0' num2str(hh)];else;hh=[num2str(hh)];end
                        if jday<10;jday=['00' num2str(jday)];elseif jday>=10 & jday<100;jday=['0' num2str(jday)];elseif jday>=100;jday=[num2str(jday)];end
                        if mi<10;mi=['0' num2str(mi)];else;mi=[num2str(mi)];end
                        if ss<10;ss=['0' num2str(ss)];else;ss=[num2str(ss)];end
                        %Create the new filename and the Output directory
                        newname=strcat(staz_ref,'.BH',component,'.',num2str(yy),'.',num2str(jday),'.',hh,mi,ss,'.mat');
                        pathout=strcat(savename,'/',num2str(jday),'/');
                        if ~exist(pathout, 'dir');mkdir(pathout);end
                        %Fill out the fields of the signal structure
                        mytrace.station=staz_ref;mytrace.channel=['BH' component];mytrace.data=z;mytrace.sampleCount=length(z);
                        mytrace.sampleRate=fs;mytrace.startTime=startTime;mytrace.endTime=endTime;mytrace.latitude=lat_coord;
                        mytrace.longitude=lon_coord;mytrace.elevation=ele_coord; 
                        %Save the seismic trace in .mat format
                        if exist([pathout ,'/', newname])==2
                            mytrace2=load([pathout ,'/', newname],'mytrace');
                            sign1=z;
                            sign2=mytrace2.mytrace.data;
                            idxzer=sign2~=0;
                            sign1(idxzer)=sign2(idxzer);
                            mytrace.data=sign1;
                            save([pathout ,'/', newname],'mytrace')
                        else
                            save([pathout ,'/', newname],'mytrace')
                        end
                end
            case 'Seed'
                %Loop through the number of files selected
                for kk=1:nsize
                    %Read and load the kk-th file .mseed
                    filename=ff{kk};
                      T = rdmseed(filename);
                        %Return the amplitude, time and info vectors
                        d=cell2mat({T(1:end).d}');
                        t=cell2mat({T(1:end).t}');
                        fs=T(1).SampleRate;nt=length(d);
                         %Conversion from counts to Volts
                        if ~isnan(bmw) & ~isempty(bmw) 
                            d = d.*bmw;
                        end
                        %Conversion from Volts to m/s
                        if ~isnan(vvel) & ~isempty(vvel) 
                            d = d./vvel;
                        end
                        %Define temporal variables and names
                        timeref=t(1);yy=year(timeref);hh=hour(timeref);mi=minute(timeref);ss=round(second(timeref));
                        timeref=datetime(timeref,'ConvertFrom','datenum');
                        timeref.Format='dd-MMM-yyyy HH:mm:SS';
                        jday=day(timeref,'dayofyear');
                        timecontrol1=datetime(datenum(yy,0,jday,hh,0,0),'ConvertFrom','datenum');
                        timecontrol1.Format='dd-MMM-yyyy HH:mm:SS';
                        timecontrol2=datetime(datenum(yy,0,jday,hh+1,0,0),'ConvertFrom','datenum');
                        timecontrol2.Format='dd-MMM-yyyy HH:mm:SS';
                        timerefend=timeref+seconds((length(d)/round(fs)));
                        timerefend.Format='dd-MMM-yyyy HH:mm:SS';
                        TimeThr=timecontrol1+seconds(1/round(fs)):seconds(1/round(fs)):timecontrol2;
                        idxtime1=find(TimeThr==timeref);
                        idxtime2=find(TimeThr==timerefend);
                        %Fill temporal gaps to create hourly files 
                        if nt~=3600*round(fs)
                            if timeref>timecontrol1
                                gapp=idxtime1;
                                d=[ones(gapp,1).*0;d];
                            end
                            if timerefend<timecontrol2
                                gapp=3600*round(fs)-idxtime2;
                                 d=[d;ones(gapp,1).*0];
                            end
                            timeref=timecontrol1;
                            mi=0;
                            ss=0;
                        end
                        startTime=datenum(timeref);
                        endTime=datenum(timeref+seconds((length(d)/round(fs)))); 
%                         endTime=t(1)+(length(d)/round(fs))/(24*3600)+(1/fs)/(24*3600);
                        if hh<10;hh=['0' num2str(hh)];else;hh=[num2str(hh)];end
                        if jday<10;jday=['00' num2str(jday)];elseif jday>=10 & jday<100;jday=['0' num2str(jday)];elseif jday>=100;jday=[num2str(jday)];end
                        if mi<10;mi=['0' num2str(mi)];else;mi=[num2str(mi)];end
                        if ss<10;ss=['0' num2str(ss)];else;ss=[num2str(ss)];end
                        %Create the new filename and the Output directory
                        newname=strcat(staz_ref,'.BH',component,'.',num2str(yy),'.',num2str(jday),'.',hh,mi,ss,'.mat');
                        pathout=strcat(savename,'/',num2str(jday),'/');
                        if ~exist(pathout, 'dir');mkdir(pathout);end
                        %Fill out the fields of the signal structure
                        mytrace.station=staz_ref;mytrace.channel=['BH' component];mytrace.data=d;mytrace.sampleCount=length(d);
                        mytrace.sampleRate=fs;mytrace.startTime=startTime;mytrace.endTime=endTime;mytrace.latitude=lat_coord;
                        mytrace.longitude=lon_coord;mytrace.elevation=ele_coord; 
                        %Save the seismic trace in .mat format
                        if exist([pathout ,'/', newname])==2
                            mytrace2=load([pathout ,'/', newname],'mytrace');
                            sign1=d;
                            sign2=mytrace2.mytrace.data;
                            idxzer=sign2~=0;
                            sign1(idxzer)=sign2(idxzer);
                            mytrace.data=sign1;
                            save([pathout ,'/', newname],'mytrace')
                        else
                            save([pathout ,'/', newname],'mytrace')
                        end
                end
            case 'Cube'
                %Control and set the Input folder
                if ~isempty(pathnam);pathname=pathname2{1};end
                %Term of research file 
                name_template = [staz_ref, '.BH',component, '.%Y.%j.%H%M%s.mseed'];
                %Create the Output directories for the conversion steps
                temp_directory=strcat(pathname,'/temp/');if ~exist(temp_directory, 'dir');mkdir(temp_directory);end
                temp2_directory=strcat(pathname,'/temp2/');if ~exist(temp2_directory, 'dir');mkdir(temp2_directory);end
                output_directory=strcat(pathname,'/convers/');if ~exist(output_directory, 'dir');mkdir(output_directory);end 
                delete(strcat(temp_directory,'*'))
                delete(strcat(temp2_directory,'*'))
                delete(strcat(output_directory,'*'))
                %Set the component/channel name  
                comps={'Z','N','E','F'};cubeformat={'*.pri0','*.pri2','*.pri1','*.pri0'};
                idxcomp=strcmp(component,comps);cubef=cubeformat{idxcomp};
                %Loop through the number of files selected
                for kk=1:nsize
                    filename=ff{kk};
                    % Covert from cube to mseed format 
                    system(['cube2mseed ', '--output-dir=' , temp_directory, ' --resample=LINEAR',' --fringe-samples=NOMINAL ', filename])
                    % Cut converted data into 1-hour or 1-day long files
                    system(['mseedcut ', '--output-dir=', temp2_directory, ' --file-length=', 'HOUR', ' ', temp_directory]);
                    % Rename files and save them in SDS-type folder/file structure.
                    system(['mseedrename ', '--template=', name_template, ' --include-pattern=', cubef, ' --transfer-mode=MOVE', ' --output-dir=', output_directory, ' ', temp2_directory] );
                    %Search all seismic traces based on the station names
                    list2=dir(strcat(output_directory,staz_ref, '.BH',component,'*.mseed'));
                    %Loop through the number of files converted
                    for jj=1:length(list2)
                        %Read and load the jj-th file .mseed
                        fileseed=strcat(output_directory,list2(jj).name);
                        T = rdmseed(fileseed);
                        %Return the amplitude, time and info vectors
                        d=cell2mat({T(1:end).d}');
                        t=cell2mat({T(1:end).t}');
                        fs=T(1).SampleRate;nt=length(d);
                        %Conversion from counts to Volts
                        if ~isnan(bmw) & ~isempty(bmw) 
                            d = d.*bmw;
                        end
                        %Conversion from Volts to m/s
                        if ~isnan(vvel) & ~isempty(vvel) 
                            d = d./vvel;
                        end
                        %Define temporal variables and names
                        timeref=t(1);yy=year(timeref);hh=hour(timeref);mi=minute(timeref);ss=round(second(timeref));
                        timeref=datetime(timeref,'ConvertFrom','datenum');
                        timeref.Format='dd-MMM-yyyy HH:mm:ss.SSSSS';
                        jday=day(timeref,'dayofyear');
                        timecontrol1=datetime(datenum(yy,0,jday,hh,0,0),'ConvertFrom','datenum');
                        timecontrol1.Format='dd-MMM-yyyy HH:mm:ss.SSSSS';
                        timecontrol2=datetime(datenum(yy,0,jday,hh+1,0,-1/round(fs)),'ConvertFrom','datenum');
                        timecontrol2.Format='dd-MMM-yyyy HH:mm:ss.SSSSS';
                        timerefend=timeref+seconds((length(d)/round(fs)))-seconds(1/round(fs));
                        timerefend.Format='dd-MMM-yyyy HH:mm:ss.SSSSS';
                        TimeThr=timecontrol1:seconds(1/round(fs)):timecontrol2;
                        idxtime1=find(TimeThr==timeref)-1;
                        idxtime2=find(TimeThr==timerefend);
                        %Fill temporal gaps to create hourly files 
                        if nt~=3600*round(fs)
                            if timeref>timecontrol1
                                gapp=idxtime1;
                                d=[ones(gapp,1).*0;d];
                            end
                            if timerefend<timecontrol2
                                gapp=3600*round(fs)-idxtime2;
                                 d=[d;ones(gapp,1).*0];
                            end
                            timeref=timecontrol1;
                            mi=0;
                            ss=0;
                        end
                        startTime=datenum(timeref);
                        endTime=datenum(timeref+seconds((length(d)/round(fs)))); 
%                         endTime=t(1)+(length(d)/round(fs))/(24*3600)+(1/fs)/(24*3600);
                        if hh<10;hh=['0' num2str(hh)];else;hh=[num2str(hh)];end
                        if jday<10;jday=['00' num2str(jday)];elseif jday>=10 & jday<100;jday=['0' num2str(jday)];elseif jday>=100;jday=[num2str(jday)];end
                        if mi<10;mi=['0' num2str(mi)];else;mi=[num2str(mi)];end
                        if ss<10;ss=['0' num2str(ss)];else;ss=[num2str(ss)];end
                        %Create the new filename and the Output directory
                        newname=strcat(staz_ref,'.BH',component,'.',num2str(yy),'.',num2str(jday),'.',hh,mi,ss,'.mat');
                        pathout=strcat(savename,'/',num2str(jday),'/');
                        if ~exist(pathout, 'dir');mkdir(pathout);end
                        %Fill out the fields of the signal structure
                        mytrace.station=staz_ref;mytrace.channel=['BH' component];mytrace.data=d;mytrace.sampleCount=length(d);
                        mytrace.sampleRate=fs;mytrace.startTime=startTime;mytrace.endTime=endTime;mytrace.latitude=lat_coord;
                        mytrace.longitude=lon_coord;mytrace.elevation=ele_coord; 
                        %Save the seismic trace in .mat format
                        if exist([pathout ,'/', newname])==2
                            mytrace2=load([pathout ,'/', newname],'mytrace');
                            sign1=d;
                            sign2=mytrace2.mytrace.data;
                            idxzer=sign2~=0;
                            sign1(idxzer)=sign2(idxzer);
                            mytrace.data=sign1;
                            save([pathout ,'/', newname],'mytrace')
                        else
                             save([pathout ,'/', newname],'mytrace')
                        end
                    end
                     %Delete the files of the previous steps
                    delete(strcat(temp_directory,'*'))
                    %Delete the files of the previous steps
                    delete(strcat(temp2_directory,'*'))
                    %Delete the files of the previous steps
                    delete(strcat(output_directory,'*'))
                    end
        end
        list3=dir(strcat(savename,'/*/',staz_ref,'.BH',component,'*.mat'));
        for gg=1:length(list3)
            filezer=strcat(list3(gg).folder,'/',list3(gg).name);
            load(filezer,'mytrace')
            d=mytrace.data;
            idx_zero=(d == 0);
            d(idx_zero) = interp1(find(~idx_zero), d(~idx_zero), find(idx_zero), 'nearest','extrap');
            mytrace.data=d;
            mytrace.units=units;
            mytrace.azimuth=az;
            mytrace.depth=depth;
            mytrace.dip=dip;
            save(filezer,'mytrace')
        end
         %Create and save the station coordinates structures
         mkcoord2(savename,sta_coord,mytrace.latitude,mytrace.longitude,mytrace.elevation)
         %Create and save the instrumental response correction structures
         mkparm2(savename,k,poles,zeros,vvel,bmw,sta_coord)
        %Command message
        set(text1, 'String', strcat('Signals',' are saved in',' OutputFolder.'));
        drawnow
        else

        %%%%%%Automatic Mode%%%%%%%%%%%%%%%%%%%%%%%
        
        savename=get(h2.savein,'string');%Get the Output directory
        pathsta=get(h2.pathsta,'string');%Get the Stations features  directory
        staz_ref=deblank(get(h2.staz_reff,'string'));
        staz_red_str=staz_ref;
        
        %Error control
        if isempty(savename);set(text1, 'String', sprintf('Error. No output directory.'));drawnow;return;end
        if isempty(pathsta);set(text1, 'String', sprintf('Error. No XML files directory.'));drawnow;return;end
        %Command message
        textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
        set(text1, 'String', textLabel); 
        drawnow
        nsize=length(ff);%Number of files selected
    

        %Load the initial signal structure
        load (strcat('mytrace.mat'))
        
        %Search all XML available
        switch cchan
            case 'Comp'
                list_xml=dir([pathsta '/*' component '.xml']);nxml=length(list_xml);
                        % Convert files on the basis of XML file availability
        for xl=1:nxml
            xml_file=[list_xml(xl).folder '/' list_xml(xl).name];
            [sta_new,~,~,~,~,~,~,~,~,~,~,~,~,~]=ReadXLMfile(xml_file);
            staz_ref=char(sta_new);            

            %Try for all tyep of components
            for ccmp=1:4
                comp_ref=comps_ref{ccmp};
                %Redefine the station name of reference
                sta_coord=staz_ref;
                 %Search specific XML file
                list_xml_new=dir([pathsta '/*'  staz_ref '*' comp_ref '.xml']);
                if ~isempty(list_xml_new)
                file_xlm_new=[list_xml_new(1).folder '/' list_xml_new(1).name];
                [~,~,lat_coord,lon_coord,ele_coord,vvel,bmw,k,az,dip,depth,units,poles,zeros]=ReadXLMfile(xml_file);
                 
         
        %Execute one of three groups of routines (.sac,.mseed,.cube files)
        switch format_file
            case 'Sac '
                %Loop through the number of files selected
                for kk=1:nsize
                    %Read and load the kk-th file .sac 
                    filename=ff{kk};
                    if ~isempty(pathnam);pathname=pathname2{kk};pathname=[pathname '/'];end
                    [~,name,exts] = fileparts(filename);
                    ix=strfind(name,staz_red_str);
                    name_cut=name(ix+length(staz_red_str):end);
                    name_cut=replace(name_cut,['H' component],['H' comp_ref]);
                    
                    file_search=[pathname '*' staz_ref '*' name_cut exts];
                    check_id=dir(file_search);
                    if ~isempty(check_id)
                    filename_new=[check_id(1).folder '/' check_id(1).name];
                    T=rsac(filename_new); 
                    %Return the amplitude, time and info vectors
                    z=T(:,2);dt=T(1,3); fs=(1/dt);yy=T(71,3);jday=T(72,3);hh=T(73,3);mi=T(74,3);ss=T(75,3)+T(76,3)/1000;nt=length(z);
                     %Conversion from counts to Volts
                        if ~isnan(bmw) & ~isempty(bmw) 
                            z = z.*bmw;
                        end
                        %Conversion from Volts to m/s
                        if ~isnan(vvel) & ~isempty(vvel) 
                            z = z./vvel;
                        end    
                    %Define temporal variables and names
                        timeref=datenum(yy,0,jday,hh,mi,ss);
                        timeref=datetime(timeref,'ConvertFrom','datenum');
%                        timeref=datetime(timeref,'ConvertFrom','datenum');
                        timeref.Format='dd-MMM-yyyy HH:mm:SS';
                        jday=day(timeref,'dayofyear');
                        timecontrol1=datetime(datenum(yy,0,jday,hh,0,0),'ConvertFrom','datenum');
                        timecontrol1.Format='dd-MMM-yyyy HH:mm:SS';
                        timecontrol2=datetime(datenum(yy,0,jday,hh+1,0,0),'ConvertFrom','datenum');
                        timecontrol2.Format='dd-MMM-yyyy HH:mm:SS';
                        timerefend=timeref+seconds((length(z)/round(fs)));
                        timerefend.Format='dd-MMM-yyyy HH:mm:SS';
                        TimeThr=timecontrol1+seconds(1/round(fs)):seconds(1/round(fs)):timecontrol2;
                        idxtime1=find(TimeThr==timeref);
                        idxtime2=find(TimeThr==timerefend);
                        %Fill temporal gaps to create hourly files 
                        if nt~=3600*round(fs)
                            if timeref>timecontrol1
                                gapp=idxtime1;
                                z=[ones(gapp,1).*0;z];
                            end
                            if timerefend<timecontrol2
                                gapp=3600*round(fs)-idxtime2;
                                 z=[z;ones(gapp,1).*0];
                            end
                            timeref=timecontrol1;
                            mi=0;
                            ss=0;
                        end
                        startTime=datenum(timeref);
                        endTime=datenum(timeref+seconds((length(z)/round(fs)))); 
%                         endTime=t(1)+(length(d)/round(fs))/(24*3600)+(1/fs)/(24*3600);
                        if hh<10;hh=['0' num2str(hh)];else;hh=[num2str(hh)];end
                        if jday<10;jday=['00' num2str(jday)];elseif jday>=10 & jday<100;jday=['0' num2str(jday)];elseif jday>=100;jday=[num2str(jday)];end
                        if mi<10;mi=['0' num2str(mi)];else;mi=[num2str(mi)];end
                        if ss<10;ss=['0' num2str(ss)];else;ss=[num2str(ss)];end
                        %Create the new filename and the Output directory
                        newname=strcat(staz_ref,'.BH',comp_ref,'.',num2str(yy),'.',num2str(jday),'.',hh,mi,ss,'.mat');
                        pathout=strcat(savename,'/',num2str(jday),'/');
                        if ~exist(pathout, 'dir');mkdir(pathout);end
                        %Fill out the fields of the signal structure
                        mytrace.station=staz_ref;mytrace.channel=['BH' comp_ref];mytrace.data=z;mytrace.sampleCount=length(z);
                        mytrace.sampleRate=fs;mytrace.startTime=startTime;mytrace.endTime=endTime;mytrace.latitude=lat_coord;
                        mytrace.longitude=lon_coord;mytrace.elevation=ele_coord; 
                        %Save the seismic trace in .mat format
                        if exist([pathout ,'/', newname])==2
                            mytrace2=load([pathout ,'/', newname],'mytrace');
                            sign1=z;
                            sign2=mytrace2.mytrace.data;
                            idxzer=sign2~=0;
                            sign1(idxzer)=sign2(idxzer);
                            mytrace.data=sign1;
                            save([pathout ,'/', newname],'mytrace')
                        else
                            save([pathout ,'/', newname],'mytrace')
                        end
                    end
                end
            case 'Seed'
                %Loop through the number of files selected
                for kk=1:nsize
                    filename=ff{kk};
                    if ~isempty(pathnam);pathname=pathname2{kk};pathname=[pathname '/'];end
                    [~,name,exts] = fileparts(filename);
                    ix=strfind(name,staz_red_str);
                    name_cut=name(ix+length(staz_red_str):end);
                    name_cut=replace(name_cut,['H' component],['H' comp_ref]);
                    
                    file_search=[pathname '*' staz_ref '*' name_cut exts];
                    check_id=dir(file_search);
                    if ~isempty(check_id)
                    filename_new=[check_id(1).folder '/' check_id(1).name];
                    %Read and load the kk-th file .mseed
                      T = rdmseed(filename_new);
                        %Return the amplitude, time and info vectors
                        d=cell2mat({T(1:end).d}');
                        t=cell2mat({T(1:end).t}');
                        fs=T(1).SampleRate;nt=length(d);
                         %Conversion from counts to Volts
                        if ~isnan(bmw) & ~isempty(bmw) 
                            d = d.*bmw;
                        end
                        %Conversion from Volts to m/s
                        if ~isnan(vvel) & ~isempty(vvel) 
                            d = d./vvel;
                        end
                        %Define temporal variables and names
                        timeref=t(1);yy=year(timeref);hh=hour(timeref);mi=minute(timeref);ss=round(second(timeref));
                        timeref=datetime(timeref,'ConvertFrom','datenum');
                        timeref.Format='dd-MMM-yyyy HH:mm:SS';
                        jday=day(timeref,'dayofyear');
                        timecontrol1=datetime(datenum(yy,0,jday,hh,0,0),'ConvertFrom','datenum');
                        timecontrol1.Format='dd-MMM-yyyy HH:mm:SS';
                        timecontrol2=datetime(datenum(yy,0,jday,hh+1,0,0),'ConvertFrom','datenum');
                        timecontrol2.Format='dd-MMM-yyyy HH:mm:SS';
                        timerefend=timeref+seconds((length(d)/round(fs)));
                        timerefend.Format='dd-MMM-yyyy HH:mm:SS';
                        TimeThr=timecontrol1+seconds(1/round(fs)):seconds(1/round(fs)):timecontrol2;
                        idxtime1=find(TimeThr==timeref);
                        idxtime2=find(TimeThr==timerefend);
                        %Fill temporal gaps to create hourly files 
                        if nt~=3600*round(fs)
                            if timeref>timecontrol1
                                gapp=idxtime1;
                                d=[ones(gapp,1).*0;d];
                            end
                            if timerefend<timecontrol2
                                gapp=3600*round(fs)-idxtime2;
                                 d=[d;ones(gapp,1).*0];
                            end
                            timeref=timecontrol1;
                            mi=0;
                            ss=0;
                        end
                        startTime=datenum(timeref);
                        endTime=datenum(timeref+seconds((length(d)/round(fs)))); 
%                         endTime=t(1)+(length(d)/round(fs))/(24*3600)+(1/fs)/(24*3600);
                        if hh<10;hh=['0' num2str(hh)];else;hh=[num2str(hh)];end
                        if jday<10;jday=['00' num2str(jday)];elseif jday>=10 & jday<100;jday=['0' num2str(jday)];elseif jday>=100;jday=[num2str(jday)];end
                        if mi<10;mi=['0' num2str(mi)];else;mi=[num2str(mi)];end
                        if ss<10;ss=['0' num2str(ss)];else;ss=[num2str(ss)];end
                        %Create the new filename and the Output directory
                        newname=strcat(staz_ref,'.BH',comp_ref,'.',num2str(yy),'.',num2str(jday),'.',hh,mi,ss,'.mat');
                        pathout=strcat(savename,'/',num2str(jday),'/');
                        if ~exist(pathout, 'dir');mkdir(pathout);end
                        %Fill out the fields of the signal structure
                        mytrace.station=staz_ref;mytrace.channel=['BH' comp_ref];mytrace.data=d;mytrace.sampleCount=length(d);
                        mytrace.sampleRate=fs;mytrace.startTime=startTime;mytrace.endTime=endTime;mytrace.latitude=lat_coord;
                        mytrace.longitude=lon_coord;mytrace.elevation=ele_coord; 
                        %Save the seismic trace in .mat format
                        if exist([pathout ,'/', newname])==2
                            mytrace2=load([pathout ,'/', newname],'mytrace');
                            sign1=d;
                            sign2=mytrace2.mytrace.data;
                            idxzer=sign2~=0;
                            sign1(idxzer)=sign2(idxzer);
                            mytrace.data=sign1;
                            save([pathout ,'/', newname],'mytrace')
                        else
                            save([pathout ,'/', newname],'mytrace')
                        end
                    end
                end
            case 'Cube'
                %Control and set the Input folder
                if ~isempty(pathnam);pathname=pathname2{1};end
                %Term of research file 
                name_template = [staz_ref, '.BH',comp_ref, '.%Y.%j.%H%M%s.mseed'];
                %Create the Output directories for the conversion steps
                temp_directory=strcat(pathname,'temp/');if ~exist(temp_directory, 'dir');mkdir(temp_directory);end
                temp2_directory=strcat(pathname,'temp2/');if ~exist(temp2_directory, 'dir');mkdir(temp2_directory);end
                output_directory=strcat(pathname,'convers/');if ~exist(output_directory, 'dir');mkdir(output_directory);end 
                delete(strcat(temp_directory,'*'))
                delete(strcat(temp2_directory,'*'))
                delete(strcat(output_directory,'*'))
                %Set the component/channel name  
                comps={'Z','N','E','F'};cubeformat={'*.pri0','*.pri2','*.pri1','*.pri0'};
                idxcomp=strcmp(component,comps);cubef=cubeformat{idxcomp};
                %Loop through the number of files selected
                for kk=1:nsize
                    filename=ff{kk};
                    [~,name,exts] = fileparts(filename);
                    ix=strfind(name,staz_red_str);
                    name_cut=name(ix+length(staz_red_str):end);
                    name_cut=replace(name_cut,['H' component],['H' comp_ref]);
                    
                    file_search=[pathname '*' staz_ref '*' name_cut exts];
                    check_id=dir(file_search);
                    if ~isempty(check_id)
                    filename_new=[check_id(1).folder '/' check_id(1).name];
                    % Covert from cube to mseed format 
                    system(['cube2mseed ', '--output-dir=' , temp_directory, ' --resample=LINEAR',' --fringe-samples=NOMINAL ', filename_new])
                    % Cut converted data into 1-hour or 1-day long files
                    system(['mseedcut ', '--output-dir=', temp2_directory, ' --file-length=', 'HOUR', ' ', temp_directory]);
                    % Rename files and save them in SDS-type folder/file structure.
                    system(['mseedrename ', '--template=', name_template, ' --include-pattern=', cubef, ' --transfer-mode=MOVE', ' --output-dir=', output_directory, ' ', temp2_directory] );
                    %Search all seismic traces based on the station names
                    list2=dir(strcat(output_directory,staz_ref, '.BH',comp_ref,'*'));
                    %Loop through the number of files converted
                    for jj=1:length(list2)
                        %Read and load the jj-th file .mseed
                        fileseed=strcat(output_directory,list2(jj).name);
                        T = rdmseed(fileseed);
                        %Return the amplitude, time and info vectors
                        d=cell2mat({T(1:end).d}');
                        t=cell2mat({T(1:end).t}');
                        fs=T(1).SampleRate;nt=length(d);
                        %Conversion from counts to Volts
                        if ~isnan(bmw) & ~isempty(bmw) 
                            d = d.*bmw;
                        end
                        %Conversion from Volts to m/s
                        if ~isnan(vvel) & ~isempty(vvel) 
                            d = d./vvel;
                        end
                        %Define temporal variables and names
                        timeref=t(1);yy=year(timeref);hh=hour(timeref);mi=minute(timeref);ss=round(second(timeref));
                        timeref=datetime(timeref,'ConvertFrom','datenum');
                        timeref.Format='dd-MMM-yyyy HH:mm:ss.SSSSS';
                        jday=day(timeref,'dayofyear');
                        timecontrol1=datetime(datenum(yy,0,jday,hh,0,0),'ConvertFrom','datenum');
                        timecontrol1.Format='dd-MMM-yyyy HH:mm:ss.SSSSS';
                        timecontrol2=datetime(datenum(yy,0,jday,hh+1,0,-1/round(fs)),'ConvertFrom','datenum');
                        timecontrol2.Format='dd-MMM-yyyy HH:mm:ss.SSSSS';
                        timerefend=timeref+seconds((length(d)/round(fs)))-seconds(1/round(fs));
                        timerefend.Format='dd-MMM-yyyy HH:mm:ss.SSSSS';
                        TimeThr=timecontrol1:seconds(1/round(fs)):timecontrol2;
                        idxtime1=find(TimeThr==timeref)-1;
                        idxtime2=find(TimeThr==timerefend);
                        %Fill temporal gaps to create hourly files 
                        if nt~=3600*round(fs)
                            if timeref>timecontrol1
                                gapp=idxtime1;
                                d=[ones(gapp,1).*0;d];
                            end
                            if timerefend<timecontrol2
                                gapp=3600*round(fs)-idxtime2;
                                 d=[d;ones(gapp,1).*0];
                            end
                            timeref=timecontrol1;
                            mi=0;
                            ss=0;
                        end
                        startTime=datenum(timeref);
                        endTime=datenum(timeref+seconds((length(d)/round(fs)))); 
%                         endTime=t(1)+(length(d)/round(fs))/(24*3600)+(1/fs)/(24*3600);
                        if hh<10;hh=['0' num2str(hh)];else;hh=[num2str(hh)];end
                        if jday<10;jday=['00' num2str(jday)];elseif jday>=10 & jday<100;jday=['0' num2str(jday)];elseif jday>=100;jday=[num2str(jday)];end
                        if mi<10;mi=['0' num2str(mi)];else;mi=[num2str(mi)];end
                        if ss<10;ss=['0' num2str(ss)];else;ss=[num2str(ss)];end
                        %Create the new filename and the Output directory
                        newname=strcat(staz_ref,'.BH',component,'.',num2str(yy),'.',num2str(jday),'.',hh,mi,ss,'.mat');
                        pathout=strcat(savename,'/',num2str(jday),'/');
                        if ~exist(pathout, 'dir');mkdir(pathout);end
                        %Fill out the fields of the signal structure
                        mytrace.station=staz_ref;mytrace.channel=['BH' component];mytrace.data=d;mytrace.sampleCount=length(d);
                        mytrace.sampleRate=fs;mytrace.startTime=startTime;mytrace.endTime=endTime;mytrace.latitude=lat_coord;
                        mytrace.longitude=lon_coord;mytrace.elevation=ele_coord; 
                        %Save the seismic trace in .mat format
                        if exist([pathout ,'/', newname])==2
                            mytrace2=load([pathout ,'/', newname],'mytrace');
                            sign1=d;
                            sign2=mytrace2.mytrace.data;
                            idxzer=sign2~=0;
                            sign1(idxzer)=sign2(idxzer);
                            mytrace.data=sign1;
                            save([pathout ,'/', newname],'mytrace')
                        else
                             save([pathout ,'/', newname],'mytrace')
                        end
                    end
                     %Delete the files of the previous steps
                    delete(strcat(temp_directory,'*'))
                    %Delete the files of the previous steps
                    delete(strcat(temp2_directory,'*'))
                    %Delete the files of the previous steps
                    delete(strcat(output_directory,'*'))
                    end
        end
        list3=dir(strcat(savename,'/*/',staz_ref,'.BH',comp_ref,'*.mat'));
        for gg=1:length(list3)
            filezer=strcat(list3(gg).folder,'/',list3(gg).name);
            load(filezer,'mytrace')
            d=mytrace.data;
            idx_zero=(d == 0);
            d(idx_zero) = interp1(find(~idx_zero), d(~idx_zero), find(idx_zero), 'nearest','extrap');
            mytrace.data=d;
            mytrace.units=units;
            mytrace.azimuth=az;
            mytrace.depth=depth;
            mytrace.dip=dip;
            save(filezer,'mytrace')
        end
        end
         %Create and save the station coordinates structures
         mkcoord2(savename,sta_coord,mytrace.latitude,mytrace.longitude,mytrace.elevation)
         %Create and save the instrumental response correction structures
         mkparm2(savename,k,poles,zeros,vvel,bmw,sta_coord)       
        end
            end
        end
         %Command message
        set(text1, 'String', strcat('Signals',' are saved in',' OutputFolder.'));
        drawnow
            case 'Chan'
                list_xml=dir([pathsta '/*.xml']);nxml=length(list_xml);
                        % Convert files on the basis of XML file availability
        for xl=1:nxml
            xml_file=[list_xml(xl).folder '/' list_xml(xl).name];
            [sta_new,comp_ref,~,~,~,~,~,~,~,~,~,~,~,~]=ReadXLMfile(xml_file);
            staz_ref=char(sta_new); 
            comp_ref=char(comp_ref);
            comp_ref=comp_ref(end);
                %Redefine the station name of reference
                sta_coord=[staz_ref '.BH' comp_ref];
                 %Search specific XML file
                list_xml_new=dir([pathsta '/*'  staz_ref '*' comp_ref '.xml']);
                if ~isempty(list_xml_new)
                file_xlm_new=[list_xml_new(1).folder '/' list_xml_new(1).name];
                [~,~,lat_coord,lon_coord,ele_coord,vvel,bmw,k,az,dip,depth,units,poles,zeros]=ReadXLMfile(xml_file);
                 
         
        %Execute one of three groups of routines (.sac,.mseed,.cube files)
        switch format_file
            case 'Sac '
                %Loop through the number of files selected
                for kk=1:nsize
                    %Read and load the kk-th file .sac 
                    filename=ff{kk};
                    if ~isempty(pathnam);pathname=pathname2{kk};pathname=[pathname '/'];end
                    [~,name,exts] = fileparts(filename);
                    ix=strfind(name,staz_red_str);
                    name_cut=name(ix+length(staz_red_str):end);
                    name_cut=replace(name_cut,['H' component],['H' comp_ref]);
                    
                    file_search=[pathname '*' staz_ref '*' name_cut exts];
                    check_id=dir(file_search);
                    if ~isempty(check_id)
                    filename_new=[check_id(1).folder '/' check_id(1).name];
                    T=rsac(filename_new); 
                    %Return the amplitude, time and info vectors
                    z=T(:,2);dt=T(1,3); fs=(1/dt);yy=T(71,3);jday=T(72,3);hh=T(73,3);mi=T(74,3);ss=T(75,3)+T(76,3)/1000;nt=length(z);
                     %Conversion from counts to Volts
                        if ~isnan(bmw) & ~isempty(bmw) 
                            z = z.*bmw;
                        end
                        %Conversion from Volts to m/s
                        if ~isnan(vvel) & ~isempty(vvel) 
                            z = z./vvel;
                        end    
                    %Define temporal variables and names
                        timeref=datenum(yy,0,jday,hh,mi,ss);
                        timeref=datetime(timeref,'ConvertFrom','datenum');
%                        timeref=datetime(timeref,'ConvertFrom','datenum');
                        timeref.Format='dd-MMM-yyyy HH:mm:SS';
                        jday=day(timeref,'dayofyear');
                        timecontrol1=datetime(datenum(yy,0,jday,hh,0,0),'ConvertFrom','datenum');
                        timecontrol1.Format='dd-MMM-yyyy HH:mm:SS';
                        timecontrol2=datetime(datenum(yy,0,jday,hh+1,0,0),'ConvertFrom','datenum');
                        timecontrol2.Format='dd-MMM-yyyy HH:mm:SS';
                        timerefend=timeref+seconds((length(z)/round(fs)));
                        timerefend.Format='dd-MMM-yyyy HH:mm:SS';
                        TimeThr=timecontrol1+seconds(1/round(fs)):seconds(1/round(fs)):timecontrol2;
                        idxtime1=find(TimeThr==timeref);
                        idxtime2=find(TimeThr==timerefend);
                        %Fill temporal gaps to create hourly files 
                        if nt~=3600*round(fs)
                            if timeref>timecontrol1
                                gapp=idxtime1;
                                z=[ones(gapp,1).*0;z];
                            end
                            if timerefend<timecontrol2
                                gapp=3600*round(fs)-idxtime2;
                                 z=[z;ones(gapp,1).*0];
                            end
                            timeref=timecontrol1;
                            mi=0;
                            ss=0;
                        end
                        startTime=datenum(timeref);
                        endTime=datenum(timeref+seconds((length(z)/round(fs)))); 
%                         endTime=t(1)+(length(d)/round(fs))/(24*3600)+(1/fs)/(24*3600);
                        if hh<10;hh=['0' num2str(hh)];else;hh=[num2str(hh)];end
                        if jday<10;jday=['00' num2str(jday)];elseif jday>=10 & jday<100;jday=['0' num2str(jday)];elseif jday>=100;jday=[num2str(jday)];end
                        if mi<10;mi=['0' num2str(mi)];else;mi=[num2str(mi)];end
                        if ss<10;ss=['0' num2str(ss)];else;ss=[num2str(ss)];end
                        %Create the new filename and the Output directory
                        newname=strcat(staz_ref,'.BH',comp_ref,'.',num2str(yy),'.',num2str(jday),'.',hh,mi,ss,'.mat');
                        pathout=strcat(savename,'/',num2str(jday),'/');
                        if ~exist(pathout, 'dir');mkdir(pathout);end
                        %Fill out the fields of the signal structure
                        mytrace.station=staz_ref;mytrace.channel=['BH' comp_ref];mytrace.data=z;mytrace.sampleCount=length(z);
                        mytrace.sampleRate=fs;mytrace.startTime=startTime;mytrace.endTime=endTime;mytrace.latitude=lat_coord;
                        mytrace.longitude=lon_coord;mytrace.elevation=ele_coord; 
                        %Save the seismic trace in .mat format
                        if exist([pathout ,'/', newname])==2
                            mytrace2=load([pathout ,'/', newname],'mytrace');
                            sign1=z;
                            sign2=mytrace2.mytrace.data;
                            idxzer=sign2~=0;
                            sign1(idxzer)=sign2(idxzer);
                            mytrace.data=sign1;
                            save([pathout ,'/', newname],'mytrace')
                        else
                            save([pathout ,'/', newname],'mytrace')
                        end
                    end
                end
            case 'Seed'
                %Loop through the number of files selected
                for kk=1:nsize
                    filename=ff{kk};
                    if ~isempty(pathnam);pathname=pathname2{kk};pathname=[pathname '/'];end
                    [~,name,exts] = fileparts(filename);
                    ix=strfind(name,staz_red_str);
                    name_cut=name(ix+length(staz_red_str):end);
                    name_cut=replace(name_cut,['H' component],['H' comp_ref]);
                    
                    file_search=[pathname '*' staz_ref '*' name_cut exts];
                    check_id=dir(file_search);
                    if ~isempty(check_id)
                    filename_new=[check_id(1).folder '/' check_id(1).name];
                    %Read and load the kk-th file .mseed
                      T = rdmseed(filename_new);
                        %Return the amplitude, time and info vectors
                        d=cell2mat({T(1:end).d}');
                        t=cell2mat({T(1:end).t}');
                        fs=T(1).SampleRate;nt=length(d);
                         %Conversion from counts to Volts
                        if ~isnan(bmw) & ~isempty(bmw) 
                            d = d.*bmw;
                        end
                        %Conversion from Volts to m/s
                        if ~isnan(vvel) & ~isempty(vvel) 
                            d = d./vvel;
                        end
                        %Define temporal variables and names
                        timeref=t(1);yy=year(timeref);hh=hour(timeref);mi=minute(timeref);ss=round(second(timeref));
                        timeref=datetime(timeref,'ConvertFrom','datenum');
                        timeref.Format='dd-MMM-yyyy HH:mm:SS';
                        jday=day(timeref,'dayofyear');
                        timecontrol1=datetime(datenum(yy,0,jday,hh,0,0),'ConvertFrom','datenum');
                        timecontrol1.Format='dd-MMM-yyyy HH:mm:SS';
                        timecontrol2=datetime(datenum(yy,0,jday,hh+1,0,0),'ConvertFrom','datenum');
                        timecontrol2.Format='dd-MMM-yyyy HH:mm:SS';
                        timerefend=timeref+seconds((length(d)/round(fs)));
                        timerefend.Format='dd-MMM-yyyy HH:mm:SS';
                        TimeThr=timecontrol1+seconds(1/round(fs)):seconds(1/round(fs)):timecontrol2;
                        idxtime1=find(TimeThr==timeref);
                        idxtime2=find(TimeThr==timerefend);
                        %Fill temporal gaps to create hourly files 
                        if nt~=3600*round(fs)
                            if timeref>timecontrol1
                                gapp=idxtime1;
                                d=[ones(gapp,1).*0;d];
                            end
                            if timerefend<timecontrol2
                                gapp=3600*round(fs)-idxtime2;
                                 d=[d;ones(gapp,1).*0];
                            end
                            timeref=timecontrol1;
                            mi=0;
                            ss=0;
                        end
                        startTime=datenum(timeref);
                        endTime=datenum(timeref+seconds((length(d)/round(fs)))); 
%                         endTime=t(1)+(length(d)/round(fs))/(24*3600)+(1/fs)/(24*3600);
                        if hh<10;hh=['0' num2str(hh)];else;hh=[num2str(hh)];end
                        if jday<10;jday=['00' num2str(jday)];elseif jday>=10 & jday<100;jday=['0' num2str(jday)];elseif jday>=100;jday=[num2str(jday)];end
                        if mi<10;mi=['0' num2str(mi)];else;mi=[num2str(mi)];end
                        if ss<10;ss=['0' num2str(ss)];else;ss=[num2str(ss)];end
                        %Create the new filename and the Output directory
                        newname=strcat(staz_ref,'.BH',comp_ref,'.',num2str(yy),'.',num2str(jday),'.',hh,mi,ss,'.mat');
                        pathout=strcat(savename,'/',num2str(jday),'/');
                        if ~exist(pathout, 'dir');mkdir(pathout);end
                        %Fill out the fields of the signal structure
                        mytrace.station=staz_ref;mytrace.channel=['BH' comp_ref];mytrace.data=d;mytrace.sampleCount=length(d);
                        mytrace.sampleRate=fs;mytrace.startTime=startTime;mytrace.endTime=endTime;mytrace.latitude=lat_coord;
                        mytrace.longitude=lon_coord;mytrace.elevation=ele_coord; 
                        %Save the seismic trace in .mat format
                        if exist([pathout ,'/', newname])==2
                            mytrace2=load([pathout ,'/', newname],'mytrace');
                            sign1=d;
                            sign2=mytrace2.mytrace.data;
                            idxzer=sign2~=0;
                            sign1(idxzer)=sign2(idxzer);
                            mytrace.data=sign1;
                            save([pathout ,'/', newname],'mytrace')
                        else
                            save([pathout ,'/', newname],'mytrace')
                        end
                    end
                end
            case 'Cube'
                %Control and set the Input folder
                if ~isempty(pathnam);pathname=pathname2{1};end
                %Term of research file 
                name_template = [staz_ref, '.BH',comp_ref, '.%Y.%j.%H%M%s.mseed'];
                %Create the Output directories for the conversion steps
                temp_directory=strcat(pathname,'temp/');if ~exist(temp_directory, 'dir');mkdir(temp_directory);end
                temp2_directory=strcat(pathname,'temp2/');if ~exist(temp2_directory, 'dir');mkdir(temp2_directory);end
                output_directory=strcat(pathname,'convers/');if ~exist(output_directory, 'dir');mkdir(output_directory);end 
                delete(strcat(temp_directory,'*'))
                delete(strcat(temp2_directory,'*'))
                delete(strcat(output_directory,'*'))
                %Set the component/channel name  
                comps={'Z','N','E','F'};cubeformat={'*.pri0','*.pri2','*.pri1','*.pri0'};
                idxcomp=strcmp(component,comps);cubef=cubeformat{idxcomp};
                %Loop through the number of files selected
                for kk=1:nsize
                    filename=ff{kk};
                    [~,name,exts] = fileparts(filename);
                    ix=strfind(name,staz_red_str);
                    name_cut=name(ix+length(staz_red_str):end);
                    name_cut=replace(name_cut,['H' component],['H' comp_ref]);
                    
                    file_search=[pathname '*' staz_ref '*' name_cut exts];
                    check_id=dir(file_search);
                    if ~isempty(check_id)
                    filename_new=[check_id(1).folder '/' check_id(1).name];
                    % Covert from cube to mseed format 
                    system(['cube2mseed ', '--output-dir=' , temp_directory, ' --resample=LINEAR',' --fringe-samples=NOMINAL ', filename_new])
                    % Cut converted data into 1-hour or 1-day long files
                    system(['mseedcut ', '--output-dir=', temp2_directory, ' --file-length=', 'HOUR', ' ', temp_directory]);
                    % Rename files and save them in SDS-type folder/file structure.
                    system(['mseedrename ', '--template=', name_template, ' --include-pattern=', cubef, ' --transfer-mode=MOVE', ' --output-dir=', output_directory, ' ', temp2_directory] );
                    %Search all seismic traces based on the station names
                    list2=dir(strcat(output_directory,staz_ref, '.BH',comp_ref,'*'));
                    %Loop through the number of files converted
                    for jj=1:length(list2)
                        %Read and load the jj-th file .mseed
                        fileseed=strcat(output_directory,list2(jj).name);
                        T = rdmseed(fileseed);
                        %Return the amplitude, time and info vectors
                        d=cell2mat({T(1:end).d}');
                        t=cell2mat({T(1:end).t}');
                        fs=T(1).SampleRate;nt=length(d);
                        %Conversion from counts to Volts
                        if ~isnan(bmw) & ~isempty(bmw) 
                            d = d.*bmw;
                        end
                        %Conversion from Volts to m/s
                        if ~isnan(vvel) & ~isempty(vvel) 
                            d = d./vvel;
                        end
                        %Define temporal variables and names
                        timeref=t(1);yy=year(timeref);hh=hour(timeref);mi=minute(timeref);ss=round(second(timeref));
                        timeref=datetime(timeref,'ConvertFrom','datenum');
                        timeref.Format='dd-MMM-yyyy HH:mm:ss.SSSSS';
                        jday=day(timeref,'dayofyear');
                        timecontrol1=datetime(datenum(yy,0,jday,hh,0,0),'ConvertFrom','datenum');
                        timecontrol1.Format='dd-MMM-yyyy HH:mm:ss.SSSSS';
                        timecontrol2=datetime(datenum(yy,0,jday,hh+1,0,-1/round(fs)),'ConvertFrom','datenum');
                        timecontrol2.Format='dd-MMM-yyyy HH:mm:ss.SSSSS';
                        timerefend=timeref+seconds((length(d)/round(fs)))-seconds(1/round(fs));
                        timerefend.Format='dd-MMM-yyyy HH:mm:ss.SSSSS';
                        TimeThr=timecontrol1:seconds(1/round(fs)):timecontrol2;
                        idxtime1=find(TimeThr==timeref)-1;
                        idxtime2=find(TimeThr==timerefend);
                        %Fill temporal gaps to create hourly files 
                        if nt~=3600*round(fs)
                            if timeref>timecontrol1
                                gapp=idxtime1;
                                d=[ones(gapp,1).*0;d];
                            end
                            if timerefend<timecontrol2
                                gapp=3600*round(fs)-idxtime2;
                                 d=[d;ones(gapp,1).*0];
                            end
                            timeref=timecontrol1;
                            mi=0;
                            ss=0;
                        end
                        startTime=datenum(timeref);
                        endTime=datenum(timeref+seconds((length(d)/round(fs)))); 
%                         endTime=t(1)+(length(d)/round(fs))/(24*3600)+(1/fs)/(24*3600);
                        if hh<10;hh=['0' num2str(hh)];else;hh=[num2str(hh)];end
                        if jday<10;jday=['00' num2str(jday)];elseif jday>=10 & jday<100;jday=['0' num2str(jday)];elseif jday>=100;jday=[num2str(jday)];end
                        if mi<10;mi=['0' num2str(mi)];else;mi=[num2str(mi)];end
                        if ss<10;ss=['0' num2str(ss)];else;ss=[num2str(ss)];end
                        %Create the new filename and the Output directory
                        newname=strcat(staz_ref,'.BH',component,'.',num2str(yy),'.',num2str(jday),'.',hh,mi,ss,'.mat');
                        pathout=strcat(savename,'/',num2str(jday),'/');
                        if ~exist(pathout, 'dir');mkdir(pathout);end
                        %Fill out the fields of the signal structure
                        mytrace.station=staz_ref;mytrace.channel=['BH' component];mytrace.data=d;mytrace.sampleCount=length(d);
                        mytrace.sampleRate=fs;mytrace.startTime=startTime;mytrace.endTime=endTime;mytrace.latitude=lat_coord;
                        mytrace.longitude=lon_coord;mytrace.elevation=ele_coord; 
                        %Save the seismic trace in .mat format
                        if exist([pathout ,'/', newname])==2
                            mytrace2=load([pathout ,'/', newname],'mytrace');
                            sign1=d;
                            sign2=mytrace2.mytrace.data;
                            idxzer=sign2~=0;
                            sign1(idxzer)=sign2(idxzer);
                            mytrace.data=sign1;
                            save([pathout ,'/', newname],'mytrace')
                        else
                             save([pathout ,'/', newname],'mytrace')
                        end
                    end
                     %Delete the files of the previous steps
                    delete(strcat(temp_directory,'*'))
                    %Delete the files of the previous steps
                    delete(strcat(temp2_directory,'*'))
                    %Delete the files of the previous steps
                    delete(strcat(output_directory,'*'))
                    end
        end
        list3=dir(strcat(savename,'/*/',staz_ref,'.BH',comp_ref,'*.mat'));
        for gg=1:length(list3)
            filezer=strcat(list3(gg).folder,'/',list3(gg).name);
            load(filezer,'mytrace')
            d=mytrace.data;
            idx_zero=(d == 0);
            d(idx_zero) = interp1(find(~idx_zero), d(~idx_zero), find(idx_zero), 'nearest','extrap');
            mytrace.data=d;
            mytrace.units=units;
            mytrace.azimuth=az;
            mytrace.depth=depth;
            mytrace.dip=dip;
            save(filezer,'mytrace')
        end
        end
         %Create and save the station coordinates structures
         mkcoord2(savename,sta_coord,mytrace.latitude,mytrace.longitude,mytrace.elevation)
         %Create and save the instrumental response correction structures
         mkparm2(savename,k,poles,zeros,vvel,bmw,sta_coord)       
                end
        end
         %Command message
        set(text1, 'String', strcat('Signals',' are saved in',' OutputFolder.'));
        drawnow
        end
        


        end
    catch
        %Error assessment
        if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end 
    end
end

%==========================================================================
%%%% Save the IRIS files in Matlab format%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save2(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    try
        t1=[];t2=[];time=[];vvel=[];bmw=[];k=[];poles=[];zeros=[];%Initialize the time and the instrumental response correction vectors 
        savename2=get(h2.irisout,'string');%Get the Output directory
        %Error control
        if isempty(savename2); set(text1, 'String', sprintf('No output directory.'));drawnow;return;end
        %Command message
        textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
        set(text1, 'String', textLabel); 
        drawnow
        %Acquisition and validation of the parameters
        sta=deblank(get(h2.sta,'string'));
        net=deblank(get(h2.net,'string'));
        loc=deblank(get(h2.loc,'string'));
        channe=deblank(get(h2.channe,'string'));
        t1=datetime(deblank(get(h2.t1,'string')));
        t2=datetime(deblank(get(h2.t2,'string')));
        time=t1:hours(1):t2+hours(1); 
        %Error control
        if isempty(t1) | isempty(t1) | isempty(time);set(text1, 'String', strcat('Error. Invalid temporal range.'));drawnow;return;end
        ntime=length(time);%Number of 1-hour seismic traces 
        %Make sure the java jar is in the path, this need only be done once per MATLAB session
        javaaddpath('IRIS-WS-2.0.19.jar'); % or later version
        %Fetch data form IRIS
        %Loop through the number of seismic traces to call back
        for kk=1:ntime-1
            %Define temporal variables and names
            tt1=datestr(time(kk),'yyyy-mm-dd HH:MM:SS' );
            tt2=datestr(time(kk+1),'yyyy-mm-dd HH:MM:SS' );
            jday=day(time(kk),'dayofyear');
            yy=num2str(year(time(kk)));hh=hour(time(kk));mi=minute(time(kk));ss=second(time(kk));
            if hh<10;hh=['0' num2str(hh)];else;hh=[num2str(hh)];end
            if jday<10;jday=['00' num2str(jday)];elseif jday>=10 & jday<100;jday=['0' num2str(jday)];elseif jday>=100;jday=[num2str(jday)];end
            if mi<10;mi=['0' num2str(mi)];else;mi=[num2str(mi)];end
            if ss<10;ss=['0' num2str(ss)];else;ss=[num2str(ss)];end 
            %Create the new filename and the Output directory
            pathout=strcat(savename2,'/',net,'/',jday,'/');if ~exist(pathout, 'dir');mkdir(pathout);end
            filemat=strcat(pathout,sta,num2str(str2num(loc)),'.BH',channe(3),'.',yy,'.',jday,'.',hh,mi,ss,'.mat') ;
            %Call back the seismic traces from IRIS database
            mytrace=irisFetch.Traces(net,sta,loc,channe,tt1,tt2,'includePZ');
            %Fill out the fields of the signal structure
            mytrace.endTime=datenum(datetime(mytrace.endTime,'ConvertFrom','datenum')+seconds(1/mytrace.sampleRate));
            mytrace.data=mytrace.data./mytrace.sensitivity;
            %Save the seismic trace in .mat format
            save(filemat,'mytrace')
        end
        %Create and save the station coordinates structures
        mkcoord(savename2,net,sta,mytrace.latitude,mytrace.longitude,mytrace.elevation,num2str(str2num(loc)))
        %Create and save the instrumental response correction structures
        sacpz=mytrace.sacpz;k=sacpz.constant;poles=sacpz.poles;zeros=sacpz.zeros;
        mkparm(savename2,net,k,poles,zeros,mytrace.sensitivity,bmw,sta,num2str(str2num(loc)))
        %Command message
        set(text1, 'String', strcat('Signal',' are saved in',' Iris OutputFolder.'));
        drawnow
    catch
       %Error assessment
       if isempty(mytrace);set(text1, 'String', strcat('Signals',' Error. Signals are not  found in',' Iris Server.'));drawnow;return;end
   end
end
end