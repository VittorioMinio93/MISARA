%%%%%%%%%%%%% Signal viewer/Instrument response correction %%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params: user parameters
% avfig: signal viewer figure
% panel_color: colors of the main panel 
% fig_color: colors of the figure
% pan1: control right side panel 
% text1: Window message
% axes1: axes of the current figure
% axis2,axes3: axes of the secondary figures
% t1_lim, t2_lim: limit of x-axis on signals axes
% slidval, slidval2: value of dynamic plots 
% dispval: displacement visualization value
% fval: filter visualization value
% h: UI control parameters
%%%%%%%%%%%h2.variables/h.variables/params.variables/variables%%%%%%%%%%%%%%%%%%%%%%%%%%
% rad,Rad: radio button within a button group 
% rad1, rad2, rad3, rad4, Rad1, Rad2: button group 
%-------------------------------- General settings
%--- General info
% sta_ref: the name of the reference station 
% component, comp: the name of the component/channel
% comp_chan: the station system used (component/channel)
% station: station names
% nsta: number of stations used
% coord: the path of the coordinate file (.mat) 
% coord_sys: the coordinate system (metrics/degrees)
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
% f1, f2: frequencies filters
% delim: type of filter
%-------------------------------- Instrumental settings
%--- Instrument correction setting
% CorrMod: Active/Deactive instrument response correction
% poles, param: the path of the instrumental response correction parameters (.mat)
% data2: the Input folder
% save2: the Output folder
% Data: the table of the instrument response correction parameters 
% stations: station names
% pols: instrument Poles from factory calibration worksheet (rad/s) 
% zers: instrument Zeroes from factory calibration worksheet (rad/s) 
% bmv : coefficient for conversion from bits to Volts (V/counts) 
% vvel: coefficient for conversion from Volts to m/s (V s/m)
% k = calibration coefficient (rad/s)
% traceout, Traceout, TRACE: uncorrected/corrected seismic trace vector/matrix
% time,temp_ref, Temporal, TEMP: time vector/matrix

function instr_corr()
%% Declaration global parameters
clc
global params
comp=params.component;
component=[];
mytrace=[];
CorrMod='Yes';
ff= cell(1);
Traceout=[]; TRACE=[];TEMP=[];temp_ref=[];stations=[];Data=[];
slidval=1;dispval=0;fval=0;fs=[];nsta=0;slidval2=1;
YLab=[];
Temporal=[];
pathname=[];pathnam=[];pathname2=[];
coord=params.coord;
param=params.poles;
axes2=[];axes3=[];
%% Signal figure
fig_color= [.94 .94 .94];
panel_color= [.97 .97 .97];

avfig= figure('Name', 'Signal Viewer', 'position',get(0,'screensize'), ...
    'NumberTitle','off','units','normalized', 'color', fig_color,'MenuBar','none');%%Create the signal figure

%Window message to user
text1= uicontrol('Style','edit','FontName','Arial','Fontsize',10,'Fontweight','bold',...
                'String','Window Message','Max',5,'Min',0,'HorizontalAlignment','left',...
                'BackgroundColor', 'w','Enable','inactive',...
                'units','normalized','Position',[.05 .02 .72 .05]);


%% Axeses
% Signal axes
axes1= axes('position',[.05 .15 .72 .75]);
hold(axes1,'on')
     grid(axes1,'on')
     grid(axes1,'minor')
     ylabel(axes1,'Amplitude')

%% Control right side panel
pan1= uipanel(avfig,'visible','on','Position',[.81 .03 .185 .94],...
    'BackgroundColor', panel_color);

uipanel(avfig,'visible','on','Position',[.81 .35 .185 .01],...
    'BackgroundColor', fig_color);

%% Text boxes
%-------------------------------- Signal setting
uicontrol(pan1,'Style','text', 'String','Staz',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .67 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','F1 (Hz)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.46 .62 .2 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','F2 (Hz)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.46 .57 .2 .02],...
    'BackgroundColor',panel_color);

%-------------------------------- Axis setting
uicontrol(pan1,'Style','text', 'String','Axis setting',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.03 .53 .4 .04],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','t1',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .50 .8 .02],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','t2',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.03 .45 .8 .02],...
    'BackgroundColor',panel_color);

%% Editable text 
%-------------------------------- Signal setting
h2.sta_ref= uicontrol(pan1,'Style','edit', 'String',params.sta_ref,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.15 .66 .3 .04],...
    'BackgroundColor','w','callback',@setlimits3);
h2.f1= uicontrol(pan1,'Style','edit', 'String','0.01',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.65 .61 .3 .04],...
    'BackgroundColor','w','callback',@setlimits2,'Enable','off');
h2.f2= uicontrol(pan1,'Style','edit', 'String','10',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.65 .56 .3 .04],...
    'BackgroundColor','w','callback',@setlimits2,'Enable','off');
h2.delim= uicontrol(pan1,'Style','popup', ...
    'String',{'Low-Pass','High-Pass','Bandpass'},'Value',1,...
    'HorizontalAlignment','left','Enable','off',...
    'Units','normalized','Position',[.65 .60 .3 .1],...
    'BackgroundColor','w','callback',@setf);

%-------------------------------- Axis setting
h2.t1_lim= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.15 .49 .3 .04],...
    'BackgroundColor','w','callback',@setlimits);
h2.t2_lim= uicontrol(pan1,'Style','edit', 'String',' ',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.15 .44 .3 .04],...
    'BackgroundColor','w','callback',@setlimits);



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

% Editable text
uicontrol(pan1,'Style','text', 'String','Correction',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.05 .75 .4 .05],...
    'BackgroundColor',panel_color);

% Instrument response correction radio buttons
h.Rad= uibuttongroup(pan1,'units','normalized','BackgroundColor',panel_color,...
    'bordertype','none','Position',[.08 .72 .9 .05]);

set(h.Rad,'SelectionChangeFcn',@rad2cbk);

h.Rad1 = uicontrol( h.Rad, 'Style','Radio','String','Yes',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.0 0 .2 1],'HandleVisibility','off');
h.Rad2 = uicontrol( h.Rad, 'Style','Radio','String','No ',...
    'BackgroundColor',panel_color, 'units','normalized','position',[.25 0 .2 1],'HandleVisibility','off');

%Set the Instrument response correction mode
if strcmp(CorrMod,'Yes')
    set (h.Rad1,'value',1);
    set (h.Rad2,'value',0);

elseif strcmp(CorrMod,'No ')
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
                'Value',0,'Position',[.815 .39 .05 .03],'callback',@logscale);
            
% Displacement visualization checkbox
uicontrol('Style','checkbox','BackgroundColor', panel_color,...
                'String','Disp/Vel','units','normalized','Tag','viscb',...
                'Value',0,'Position',[.865 .39 .05 .03],'callback',@dispvisuale);
            
% Filter visualization checkbox
uicontrol('Style','checkbox','BackgroundColor', panel_color,...
                'String','Filter','units','normalized','Tag','filt',...
                'Value',0,'Position',[.835 .58 .04 .03],'callback',@useFilter);

% Editable texts
uicontrol(pan1,'Style','text', 'String','Signal Viewer/Signal correction',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.03 .24 .8 .07],...
    'BackgroundColor',panel_color);

% Calculation and save button
uicontrol('style','pushbutton', 'string','Plot signals','units','normalized',...
    'position',[.815 .22 .08 .0494], 'callback',@calc); 
uicontrol('style','pushbutton', 'string','Save Plot','units','normalized',...
    'position',[.815 .14 .08 .0494], 'callback',@save1);
uicontrol('style','pushbutton', 'string','Plot all stations','units','normalized',...
    'position',[.905 .22 .08 .0494], 'callback',@save2); 
uicontrol('style','pushbutton', 'string','Plot Response','units','normalized',...
    'position',[.905 .14 .08 .0494], 'callback',@save3); 
uicontrol('style','pushbutton', 'string','Save Data','units','normalized',...
    'position',[.815 .06 .08 .0494], 'callback',@savecorr);

%% Other functions (nested)
%%%%Return the list of files manually selected%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function openfnc(~,~)
    Traceout=[]; TRACE=[];TEMP=[];temp_ref=[];stations=[];Data=[];%Initialize the signals vectors
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);%Clear the signal axes
    sta_ref=get(h2.sta_ref,'string');%Get the station name of reference
    tt= 0;pathnam=[];pathname2=[];
    ff=cell(1);%Initialize the list of files
    search_sac=strcat(sta_ref,'.BH',comp,'*.mat');
     component=comp;
    [fname, pathname] = uigetfile({search_sac;'*.*'},'File Selector',...
        params.data2,'MultiSelect','on');%Selection of the files     
    
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
     Traceout=[]; TRACE=[];TEMP=[];temp_ref=[];stations=[];Data=[];%Initialize the signals vectors
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1);%Clear the signal axes
    pathname=[];
    pathnam= uigetdir(params.data2);
    pathnam=[pathnam '/'];
    sta_ref=deblank(get(h2.sta_ref,'string'));%Get the station name of reference
    search_sac=strcat(sta_ref,'.BH',comp,'*.mat');
    component=comp;
    
    if pathnam~=0
    % Get all files in the directory and subdirectories
    [ff,pathname2]= getAllFiles(pathnam,search_sac);
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
    
end

%==========================================================================
%%%%Return the instrument response mode%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rad2cbk(~,eventdata)
    
    CorrMod = get(eventdata.NewValue,'String');   
    
end

%==========================================================================
%%%%Plot the main seismic trace%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calc(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1)%Clear the axes
    try
        sta_ref=deblank(get(h2.sta_ref,'string'));%Get the station name of reference
        %Validation of the frequencies filters
        if fval==1
            f1=str2num(get(h2.f1,'string'));if f1<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
            f2=str2num(get(h2.f2,'string'));if f2<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
            delim=get(h2.delim,'Value');if delim==3 & f2<=f1;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        end
        %Command message
        textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
        set(text1, 'String', textLabel); 
        drawnow
        %Reading .mat files
        [~, len]=size (ff);%Number of files selected
        Traceout=[];%Initialize the signal vector
        Temporal=[];%Initialize the time vector
        Data=[];%Initialize the table of the instrument response correction parameters
        %Active/Deactive the instrument response correction mode
        switch CorrMod
            case 'Yes' 
                load(param)
                stations=Data.stations;

                %Redefine the station name of reference
                switch params.comp_chan
                    case 'Comp'
                        sta=sta_ref;   
                    case 'Chan'
                        sta=[sta_ref '.BH' comp];
                end
                %Check instrument response parameters
                idx_sta=strcmp(stations,sta);
                sta_check=stations(idx_sta);
                if isempty(sta_check);set(text1, 'String', 'Error. No parameters file directory or No correct parameters file.');return;end
                %Select the analysis parameters
                idx=strcmp(stations,sta);
                pols=Data.pols{idx};
                zers=Data.zers{idx};
                bmw=Data.bmw(idx);
                k=Data.k(idx);
                vvel=Data.vvel(idx);
                %Loop through the number of files selected
                for i=1:len
                    %Load the i-th file
                    file_name=ff{i};
                    load(file_name); 
                    z=mytrace.data;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate; nyq=0.5/dt; nt=mytrace.sampleCount;
                    startTime=mytrace.startTime;endTime=mytrace.endTime;
                    %Rescale and remove the best straight-fit line
                    z=z-mean(z);
                    z=detrend(z);
                    %Create the time vector
                    Time1=datetime(startTime,'ConvertFrom','datenum');%Convert in datetime format
                    Time2=datetime(endTime,'ConvertFrom','datenum');%Convert in datetime format
                    Time1.Format = 'dd-MM-yyyy HH:mm:ss.SSSSS';
                    Time2.Format = 'dd-MM-yyyy HH:mm:ss.SSSSS';
                    time=Time1+seconds(dt):seconds(dt):Time2;
                    time=datenum(time');
                    %Apply the instrument response correction
                    traceout = instrument_resp_remove(nyq,pols,zers,k,bmw,vvel,z);
                    %Eventually, filter the corrected seismic traces
                    if fval==1;traceout=settypeFilter2(traceout,delim,f1,f2,fs);end
                    %Concatenate the results
                    Traceout=[Traceout;traceout];
                    Temporal=[Temporal;time];            
                end
                %Set the label of the y-axis
                YLab='Amplitude corr.';
                %Plot the corrected seismic traces
                ploting1(axes1)
            case 'No '  
                %Loop through the number of files selected
                for i=1:len
                    %Load the i-th file
                    file_name=ff{i};
                    load(file_name)
                    z=mytrace.data;dt=1/mytrace.sampleRate; fs=mytrace.sampleRate; nyq=0.5/dt; nt=mytrace.sampleCount;
                    startTime=mytrace.startTime;endTime=mytrace.endTime;
                    %Rescale and remove the best straight-fit line
                    z=z-mean(z);
                    z=detrend(z);
                    %Create the time vector
                    Time1=datetime(startTime,'ConvertFrom','datenum');%Convert in datetime format
                    Time2=datetime(endTime,'ConvertFrom','datenum');%Convert in datetime format
                    Time1.Format = 'dd-MM-yyyy HH:mm:ss.SSSSS';
                    Time2.Format = 'dd-MM-yyyy HH:mm:ss.SSSSS';
                    time=Time1+seconds(dt):seconds(dt):Time2;
                    time=datenum(time');
                    traceout =z;
                    %Eventually, filter the corrected seismic traces
                    if fval==1;traceout=settypeFilter2(traceout,delim,f1,f2,fs);end
                    %Concatenate the results
                    Traceout=[Traceout;traceout];
                    Temporal=[Temporal;time];           
                end
                %Set the label of the y-axis
                YLab='Amplitude';
                %Plot the seismic traces
                ploting1(axes1)       
        end
    catch
        %Error assessment
        if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
        if isempty(Data);set(text1, 'String','Error. No parameters file directory or No correct parameters file.');drawnow;return;end
        if isempty(z);set(text1, 'String','Error. No correct signals data.');drawnow;return;end  
        if f1>fs/2 | f2>fs/2;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    end
    %Command message
    set(text1, 'String','Calculation finished.');
    drawnow
end

%==========================================================================
%%%%Plot of the seismic trace%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting1(ax)
     try
         %Redefine the station name of reference
         switch params.comp_chan 
             case'Comp'
                 sta_ref=get(h2.sta_ref,'string');
             case 'Chan'
                 sta_ref=strcat(get(h2.sta_ref,'string'),'.BH',comp);   
         end
         cla(ax)%Clear axes
         station=sta_ref;
         %Plot the signal in velocity/displacement
         if dispval==0
             traceout=Traceout;
         elseif dispval==1
             traceout=cumtrapz(Traceout);    
         end
         time=datetime(Temporal,'ConvertFrom','datenum');%Convert in datetime format
         %Plot the normalized amplitude of the seismic trace
         plot(ax,time,traceout./max(traceout))
         legend(ax,station)
         ylabel(ax,YLab)
         ylim(ax,[-1 1])
         xlim(ax,[min(time) max(time)])
     catch
         %Error assessment
         if size(Traceout,1)~=size(Temporal,1);set(text1, 'String','Error. Incoherent dimensions of the Signals array.');drawnow;return;end
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

    else
        %Command message
        set(text1, 'String', 'Linear axis chosen.');
        %Set the linear scale
        set(axes1,'yScale','linear');
    end
end

%==========================================================================
%%%%Set the limits of the x-axis on the coheregram and rms axes%%%%%%%%%%%%
function setlimits(~,~)
    try
        %Manual settings
        %Acquisition time limits
        t1_lim=datetime(get(h2.t1_lim,'string'),'InputFormat','yyMMdd-HHmmSS');%Convert in datetime format
        t2_lim=datetime(get(h2.t2_lim,'string'),'InputFormat','yyMMdd-HHmmSS'); %Convert in datetime format
        xlim(axes1,[t1_lim t2_lim])
    catch
        %Default settings
        if ~isempty(cell2mat(Temporal)) 
            time=datetime(Temporal,'ConvertFrom','datenum');%Convert in datetime format
            xlim(axes1,[min(time) max(time)])
        end       
    end
end

%==========================================================================
%%%%Control and reset automatically frequencies filters%%%%%%%%%%%%%%%%%%%%
function setlimits2(~,~)
    %Acquisition frequencies
    f1=str2num(get(h2.f1,'string'));
    f2=str2num(get(h2.f2,'string'));
    
    %Control and eventually reset frequencies
    if (size(f1,1)==0);set(h2.f1,'string','0.01');end
    if (size(f2,1)==0);set(h2.f2,'string','10');end
end

%==========================================================================
%%%%Control and reset automatically frequencies filters%%%%%%%%%%%%%%%%%%%%
function setlimits3(~,~)
     %Acquisition station name of reference
     sta_ref=deblank(get(h2.sta_ref,'string'));
     ff=cell(1);%Initialize the list of files
     
     %Control and eventually reset station name of reference
     if (size(sta_ref,2)==0);set(h2.sta_ref,'string',params.sta_ref);end
end

%==========================================================================
%%%%%%%%Save Signal plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save1(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Acquisition and validation of the parameters
    f1=str2num(get(h2.f1,'string'));f2=str2num(get(h2.f2,'string'));
    if isempty(Temporal) | isempty (Traceout);set(text1, 'String','Error. Signals matrix is empty.');drawnow;return;end
    if f1>fs/2 | f2>fs/2;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    try
        sta_ref=deblank(get(h2.sta_ref,'string'));%Get the station name of reference
        %Redefine the station name of reference
        switch params.comp_chan 
            case'Comp'
                station=sta_ref;
            case 'Chan'
                station=strcat(sta_ref,'.BH',comp);  
        end
        %Normalize the seismic traces
        traceout=Traceout;
        traceout=traceout./max(traceout);
        time=datetime(Temporal,'ConvertFrom','datenum');%Convert in datetime format
        %Plot the normalized amplitude of the seismic trace
        new_figure=figure('Name','Signal plot','NumberTitle','off');
        plot(time,traceout)
        legend(station)
        ylabel(YLab)
        ylim([-1 1])
        xlim([min(time) max(time)])
        grid on; grid minor;
    catch
        %Error assessment
        if size(Temporal,1)~=size(Traceout,1);set(text1, 'String','Error. Incoherent dimensions of the Signals array.');drawnow;return;end
    end    
end

%==========================================================================
%%%%%%%%Find all seismic traces available%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save2(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    try
        %Validation of the frequencies filters
        if fval==1
            f1=str2num(get(h2.f1,'string'));if f1<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
            f2=str2num(get(h2.f2,'string'));if f2<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
            delim=get(h2.delim,'Value');if delim==3 & f2<=f1;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        end
        %Command message
        textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
        set(text1, 'String', textLabel); 
        drawnow
        %Redefine the station name of reference
        switch params.comp_chan 
            case 'Comp'
                sta_ref=deblank(get(h2.sta_ref,'string'));
            case 'Chan'
                sta_ref=strcat(get(h2.sta_ref,'string'),'.BH',comp);   
        end
        %Active/Deactive the instrument response correction mode
        switch CorrMod
            case 'Yes'
                %Reading .mat files
                [~, len]=size (ff);%Number of files selected
                TRACE=[];%Initialize the signals matrix
                TEMP=[];%Initialize the time matrix
                %Loop through the number of files selected
                for kk=1:len
                    filez=ff{kk};%Return the kk-th files
                    %Control and set the Input folder
                    if ~isempty(pathnam);pathname=pathname2{kk};pathname=[pathname '/'];end
                    %Find and correct all seismic traces available
                    [TraceOut,TimeOut,nsta,stations]=plotAllSta(pathname,filez,param,component,params.comp_chan,params.coord,text1);
                    %Eventually, filter the corrected seismic traces
                    if fval==1;TraceOut=settypeFilter(TraceOut,delim,f1,f2,fs);end
                    idx=strcmp(stations,sta_ref);
                    %Concatenate the results
                    TRACE=[TRACE;TraceOut];
                    TEMP=[TEMP;TimeOut];
                end
                temp_ref=datetime(cell2mat(TEMP(:,idx)),'ConvertFrom','datenum');%Convert in datetime format                
                %Prepare the new figure
                new_figure=figure('units','normalized','outerposition',[0 0 1 1],'Name','Plot all stations','NumberTitle','off');
                axes2= axes('position',[.08 .12 .83 .82]);
                hold(axes2,'on')
                grid(axes2,'on')
                grid(axes2,'minor')
                %Next and pre button
                uicontrol('style','pushbutton', 'string','Next','units','normalized',...
                    'position',[.925 .6 .06 .04], 'callback',@sign2);
                uicontrol('style','pushbutton', 'string','Pre','units','normalized',...
                    'position',[.925 .66 .06 .04], 'callback',@sign1);
                %Set the block of signals
                slidval=nsta;
                %Plot all seismic corrected traces available
                ploting2(slidval,axes2)
                %Command message
                set(text1, 'String','Calculation finished.');
                drawnow
            case 'No ' 
                %Reading .mat files
                [~, len]=size (ff);%Number of files selected
                TRACE=[];%Initialize the signals matrix
                TEMP=[];%Initialize the time matrix
                %Loop through the number of files selected
                for kk=1:len
                    filez=ff{kk};%Return the kk-th files
                    %Control and set the Input folder
                    if ~isempty(pathnam);pathname=pathname2{kk};pathname=[pathname '/'];end
                     %Find all seismic traces available
                    [TraceOut,TimeOut,nsta,stations]=plotAllStaNoCorr(pathname,filez,coord,component,params.comp_chan,text1);
                    %Eventually, filter the seismic traces
                    if fval==1;TraceOut=settypeFilter(TraceOut,delim,f1,f2,fs);end
                    idx=strcmp(stations,sta_ref);
                    %Concatenate the results
                    TRACE=[TRACE;TraceOut];
                    TEMP=[TEMP;TimeOut];
                end
                temp_ref=datetime(cell2mat(TEMP(:,idx)),'ConvertFrom','datenum');%Convert in datetime format                               
                %Prepare the new figure
                new_figure=figure('units','normalized','outerposition',[0 0 1 1],'Name','Plot all stations','NumberTitle','off');
                axes2= axes('position',[.08 .12 .83 .82]);
                hold(axes2,'on')
                grid(axes2,'on')
                grid(axes2,'minor')
                %Next and pre button
                uicontrol('style','pushbutton', 'string','Next','units','normalized',...
                    'position',[.925 .6 .06 .04], 'callback',@sign2);
                uicontrol('style','pushbutton', 'string','Pre','units','normalized',...
                    'position',[.925 .66 .06 .04], 'callback',@sign1);
                %Set the block of signals
                slidval=nsta;
                %Plot all seismic traces available
                ploting2(slidval,axes2) 
                %Command message
                set(text1, 'String','Calculation finished.');
                drawnow
        end
    catch
        %Error assessment
        if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
        if isempty(params.coord);set(text1, 'String','Error. No coordinate file directory or No correct coordinate file.');drawnow;return;end
        if f1>fs/2 | f2>fs/2;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
        if isempty(Data);set(text1, 'String','Error. No parameters file directory or No correct parameters file.');drawnow;return;end
    end
end

%==========================================================================
%%%%%%%% Recall The Bode plot function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save3(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Error control
    if isempty(param);set(text1, 'String','Error. No parameters file directory or No correct parameters file.');drawnow;return;end
    %Prepare the new figure and action button
    new_figure=figure('Name','Bode plot','NumberTitle','off');
    axes3=axes;
    %Next and pre button
    uicontrol('style','pushbutton', 'string','Next','units','normalized',...
    'position',[.925 .6 .06 .04], 'callback',@ssign2);
    uicontrol('style','pushbutton', 'string','Pre','units','normalized',...
    'position',[.925 .66 .06 .04], 'callback',@ssign1);
   %Show the Bode plot
    ploting3(slidval2,axes3)
end

%==========================================================================
%%%%%%%% Plot all seismic traces available%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting2(ii,ax)
    try
        cla(ax)%Clear axes       
        trace=cell2mat(TRACE(:,ii));%Select the ii-th block of the signals matrix
        %Plot the signal in velocity/displacement
        if dispval==1
            trace=cumtrapz(trace);
        end
        temp=datetime(cell2mat(TEMP(:,ii)),'ConvertFrom','datenum');%Convert the ii-th block of the time matrix in datetime format
        %Plot the normalized amplitude of the seismic trace
        plot(ax,temp,trace./max(trace))
        title(ax,stations(ii))
        % x-axis limits control
        if size(temp,1)~=0;xlim([min(temp_ref) max(temp_ref)]);end
    catch
        %Error assessment
        if size(trace,1)~=size(temp,1);set(text1, 'String','Error. Incoherent dimensions of the Signals array.');drawnow;return;end
    end    
end

%==========================================================================
%%%%%%%%%%%%%%%%%%%%The Bode plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting3(ii,ax)
    %Load the instrument response correction parameters
    load(param)
    %Select the analysis parameters
    stations=Data.stations;nsta=length(stations);
    idx=ii;
    pols=Data.pols{idx};
    zers=Data.zers{idx};
    k=Data.k(idx);
    %Create the dynamical system model
    sys=zpk(zers,pols,k);
    %Plot the frequency response
    bode(ax,sys)
    grid (ax,'on')
    grid(ax,'minor')
    title(stations{idx})
end

%==========================================================================
%%%%Refresh the Signals axes through the pre button%%%%%%%%%%%%%%%%%%%%%%%%
function sign1(~,~)
   %Error control
   if isempty(TRACE);set(text1, 'String','Error. Signals matrix is empty.');drawnow;return;end
   %Decrease the slidval value
   slidval=slidval-1;if slidval<=0;slidval=1;end
   %Plot the new block of signals matrix 
   ploting2(slidval,axes2);
        
end

%==========================================================================
%%%%Refresh the Signals axes through the next button%%%%%%%%%%%%%%%%%%%%%%%
function sign2(~,~)
   %Error control
   if isempty(TRACE);set(text1, 'String','Error. Signals matrix is empty.');drawnow;return;end
   %Increase the slidval value
   slidval=slidval+1;if slidval>nsta;slidval=nsta;end
   %Plot the new block of signals matrix 
   ploting2(slidval,axes2);
        
end

%==========================================================================
%%%%Refresh the parameter selection through the pre button%%%%%%%%%%%%%%%%%
function ssign1(~,~)
   %Decrease the slidval2 value
   slidval2=slidval2-1;if slidval2<=0;slidval2=1;end
   %Show the Bode plot of the new block of parameters
   ploting3(slidval2,axes3)
        
end

%==========================================================================
%%%%Refresh the parameter selection through the next button%%%%%%%%%%%%%%%%
function ssign2(~,~)
   %Increase the slidval2 value
   slidval2=slidval2+1;if slidval2>nsta;slidval2=nsta;end 
   %Show the Bode plot of the new block of parameters
   ploting3(slidval2,axes3)
end

%==========================================================================
%%%%%%%%%%%%%%%%%%%%%Active/Deactive filter mode%%%%%%%%%%%%%%%%%%%%%%%%%%%
function useFilter(h, ~)
    if get(h,'Value')==get(h,'Max')
        if get(h2.delim,'Value')==1
            %Low-pass filter
            set(h2.delim,'Enable','on');
            set(h2.f1,'Enable','off');
            set(h2.f2,'Enable','on');
        elseif get(h2.delim,'Value')==2
            %High-pass filter 
            set(h2.delim,'Enable','off');
            set(h2.f1,'Enable','on');
            set(h2.f2,'Enable','off');
          elseif get(h2.delim,'Value')==3
            %Band-pass filter
            set(h2.delim,'Enable','on');
            set(h2.f1,'Enable','on');
            set(h2.f2,'Enable','on');
        end
        %Activation filter mode
        fval=1;
    else
         %No filter
         set(h2.delim,'Enable','off');
         set(h2.f1,'Enable','off');
         set(h2.f2,'Enable','off');
         %Dectivation filter mode
         fval=0;
    end   
end

%==========================================================================
%%%%%%%%%%%%%%%%%%%%%Set the frequencies filters %%%%%%%%%%%%%%%%%%%%%%%%%%
function setf(h, ~)
    if get(h,'Value')==1
        %Low-pass filter
        set(h2.f1,'Enable','off');
        set(h2.f2,'Enable','on');
    elseif get(h,'Value')==2
        %High-pass filter
        set(h2.f1,'Enable','on');
        set(h2.f2,'Enable','off');
    elseif get(h2.delim,'Value')==3
        %Band-pass filter
        set(h2.f1,'Enable','on');
        set(h2.f2,'Enable','on');      
    end  
end

%==========================================================================
%%%%%%%%%%%%%%%%%%%%%Set the visualization mode %%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispvisuale(h, ~)
    if ~isempty(Traceout)
    if (get(h,'Value') == get(h,'Max'))
        %Command message
        set(text1, 'String', 'Plotting in displacement.');
        %Set the visualization mode value
        dispval=1;
        % Plot the seismic traces in displacement
        ploting1(axes1)
    else
        %Command message
        set(text1, 'String', 'Plotting in velocity.');
        %Set the visualization mode value
        dispval=0;
        % Plot the seismic traces in velocity
        ploting1(axes1)
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
%%%%%%%%%%%%%%%%%%%%%Save the corrected seismic traces %%%%%%%%%%%%%%%%%%%%
function savecorr(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Validation of the directories
    if isempty(params.save2);set(text1, 'String', sprintf('Error. No output directory.'));drawnow;return;end  
    if isempty(cell2mat(ff));set(text1, 'String','Error. No files are loaded.');drawnow;return;end
    %Validation of the results
    if (~isempty(param)) & (CorrMod=='Yes')
        %Reading .mat files
        [~, len]=size (ff);%Number of files selected
        %Loop through the number of files selected
        for kk=1:len
            filez=ff{kk};%Return the kk-th file
            %Control and set the Input folder
            if ~isempty(pathnam);pathname=pathname2{kk};pathname=[pathname '/'];end
            %Save the kk-th corrected seismic trace
            savenewfile(pathname,filez,param,params.save2,text1,params.comp_chan,component,params.coord)
        end
    elseif (~isempty(param)) & (CorrMod=='No ')
        %Error message
        set(text1, 'String', sprintf('Error. Correction mode is no activated.'))
        drawnow
        return
    elseif (isempty(param))
        %Error message
        set(text1, 'String','Error. No parameters file directory or No correct parameters file.');
        drawnow
        return
    end
end
end