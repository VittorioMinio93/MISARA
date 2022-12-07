%%%%%%%%%%%%% Viewing of a location map and Beam Pattern analysis %%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params: user parameters
% Fig1: Map and Beam Pattern figure
% panel_color: colors of the main panel 
% fig_color: colors of the figure
% pan1: control right side panel 
% text1: Window message
% axes1,axes2: axes of the current figure
% slidval: value of dynamic plots 
% h: UI control parameters
%%%%%%%%%%%h.variables/params.variables/variables%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------- General settings
%--- General info
% station: station names
% map: the path of the DEM file (.tif)
% coord: the path of the coordinate file (.mat) 
%-------------------------------- Instrumental settings
%--- Coordinate system/Beam Pattern setting
% coord_sys: the coordinate system (metrics/degrees)
% fmin: the minimum frequency for the Beam pattern analysis (Hz)
% fmax: the maximum frequency for the Beam pattern analysis(Hz)
% step: the frequency step for the Beam pattern analysis (Hz)
% F: range of frequency analysis (Hz)  
% smin: the minimum slowness for the Beam pattern analysis (s/km)
% smax: the miximum slowness for the Beam pattern analysis (s/km)
% ns: the dimensions of the square grid slowness for the Beam pattern analysis
% sx/Sx,sy,Sy: the slowness grid vectors
% B,BeamP: the beam pattern matrix
% zmin, zmax, zstep: the elevation limits and step (km)
% data,dataCoord: the table of the station coordinates
% station: station names
% xs, ys, zs, utm_zone: coordinates in UTM system
% LAT, LON, ELE: coordinates in wgs84 system
% X,Y,A: georeferenced DEM grid NxM (UTM system/WGS84 system)
% xLimits, yLimits: georeferenced limits of the DEM grid in m/decimal degrees

function arraymap()
%% Declaration global parameters
clc
global params dataCoord
station=[];utm_zone=[];
A=[];
X=[];
Y=[];
xx=[];
yy=[];
xs=[];
ys=[];
zs=[];
LAT=[];
LON=[];
ELE=[];
BeamP={};
Sx=[];
Sy=[];
xLimits=[];
yLimits=[];
slidval=1;
F=[];smin=[];smax=[];ns=[];fmin=[];fmax=[];step=[];
zmin=[];zmax=[];zstep=[];utm_zone=[];

%% Map and Beam Pattern figure
fig_color= [.94 .94 .94];
panel_color= [.97 .97 .97];


Fig1= figure('Name', 'Beam Pattern', 'position',get(0,'screensize'), ...
    'NumberTitle','off','units','normalized', 'color', fig_color,'MenuBar','none');%%Create the map and beam pattern figure

%Window message to user
text1= uicontrol('Style','edit','FontName','Arial','Fontsize',10,'Fontweight','bold',...
                'String','Window Message','Max',5,'Min',0,'HorizontalAlignment','left',...
                'BackgroundColor', 'w','Enable','inactive',...
                'units','normalized','Position',[.05 .015 .3 .16]);
%% Axeses
% Map axes
axes1= axes('position',[.16 .28 .35 .7]);
hold(axes1,'on')
colorbar(axes1)
colormap(axes1,'parula')
%Select the correct axis labels
switch params.coord_sys
    case 'm'
   xlabel(axes1,'X (m)')
   ylabel(axes1,'Y (m)')
    case 'degrees'
    ylabel(axes1,'Latitude (decimal degrees)')
    xlabel(axes1,'Longitude (decimal degrees)')       
end
axis(axes1,'equal');
axis(axes1,'tight');

% Beam pattern axes
axes2= axes('position',[0.71 0.50 0.28 0.47]);
hold(axes2,'on')
colormap(axes2,'parula')
colorbar(axes2)
xlabel(axes2,'S_x (s/km)', 'fontsize', 8);
ylabel(axes2,'S_y (s/km)', 'fontsize', 8);

%% Control right side panel
pan1= uipanel(Fig1,'visible','on','Position',[.67 .01 .33 .42],...
    'BackgroundColor', panel_color);
uipanel(pan1,'visible','on','Position',[.01 .48 .98 .01],...
    'BackgroundColor', fig_color);

%% Text boxes
%-------------------------------- Map setting
uicontrol(pan1,'Style','text', 'String','Station Map setting',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.7 .88 .4 .07],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Zmin (km)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.7 .78 .2 .07],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Zmax (km)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.7 .665 .2 .07],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Step (km)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.7 .55 .2 .07],...
    'BackgroundColor',panel_color);

%-------------------------------- Beam Pattern setting
uicontrol(pan1,'Style','text', 'String','Beam Pattern setting',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.1 .88 .4 .07],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Fmin (Hz)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.1 .78 .2 .07],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Fmax (Hz)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.1 .665 .2 .07],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Step (s/km)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.1 .55 .2 .07],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Smin (s/km)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.38 .78 .2 .07],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Smax (s/km)',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.38 .665 .2 .07],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Npoint',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.38 .55 .2 .07],...
    'BackgroundColor',panel_color);

%% Editable text 
%-------------------------------- Map setting
h.zmin= uicontrol(pan1,'Style','edit', 'String','1.0',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.84 .79 .1 .07],...
    'BackgroundColor','w','callback',@deflimits);
h.zmax= uicontrol(pan1,'Style','edit', 'String','3.3',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.84 .67 .1 .07],...
    'BackgroundColor','w','callback',@deflimits);
h.zstep= uicontrol(pan1,'Style','edit', 'String','0.05',...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.84 .55 .1 .07],...
    'BackgroundColor','w','callback',@deflimits);

%-------------------------------- Beam Pattern setting
h.fmin= uicontrol(pan1,'Style','edit', 'String',params.fmin,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.24 .79 .1 .07],...
    'BackgroundColor','w','callback',@deflimits);
h.fmax= uicontrol(pan1,'Style','edit', 'String',params.fmax,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.24 .67 .1 .07],...
    'BackgroundColor','w','callback',@deflimits);
h.step= uicontrol(pan1,'Style','edit', 'String',params.step,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.24 .55 .1 .07],...
    'BackgroundColor','w','callback',@deflimits);
h.smin= uicontrol(pan1,'Style','edit', 'String',params.smin,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.54 .79 .1 .07],...
    'BackgroundColor','w','callback',@deflimits);
h.smax= uicontrol(pan1,'Style','edit', 'String',params.smax,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.54 .67 .1 .07],...
    'BackgroundColor','w','callback',@deflimits);
h.ns= uicontrol(pan1,'Style','edit', 'String',params.ns,...
    'HorizontalAlignment','left',...
    'Units','normalized','Position',[.54 .55 .1 .07],...
    'BackgroundColor','w','callback',@deflimits);

%% Buttons
% Editable texts
uicontrol(pan1,'Style','text', 'String','Beam Pattern',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.1 .35 .2 .07],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Stations map',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.7 .35 .2 .07],...
    'BackgroundColor',panel_color);
uicontrol(pan1,'Style','text', 'String','Coordinates',...
    'HorizontalAlignment','left','Fontweight','bold',...
    'Units','normalized','Position',[.4 .35 .3 .07],...
    'BackgroundColor',panel_color);

% Calculation, upload and save button
uicontrol('style','pushbutton', 'string','Calculation','units','normalized',...
    'position',[.700 .1 .08 .0494], 'callback',@calc);
uicontrol('style','pushbutton', 'string','Calculation','units','normalized',...
    'position',[.898 .1 .08 .0494], 'callback',@calc2);
uicontrol('style','pushbutton', 'string','Save','units','normalized',...
    'position',[.700 .03 .08 .0494], 'callback',@save1);
uicontrol('style','pushbutton', 'string','Save','units','normalized',...
    'position',[.898 .03 .08 .0494], 'callback',@save2);
uicontrol('style','pushbutton', 'string','Upload','units','normalized',...
    'position',[.799 .1 .08 .0494], 'callback',@upCoord);
uicontrol('style','pushbutton', 'string','Save','units','normalized',...
    'position',[.799 .03 .08 .0494], 'callback',@savCoord);

% Next and pre button
uicontrol('style','pushbutton', 'string','Next','units','normalized',...
    'position',[.66 .85 .03 .04], 'callback',@beam2);
uicontrol('style','pushbutton', 'string','Pre','units','normalized',...
    'position',[.66 .9 .03 .04], 'callback',@beam1);

%% Tables
% Table of stations coordinates
try
 coord_sys=params.coord_sys;load(params.coord)%Load coordinates of the stations
 columnformat = {'logical','numeric', 'numeric','numeric','char'};%Set the columns format of the table
    
    %Load the table on the basis of the geographic system (wgs84/Utm )
    switch coord_sys
        case 'degrees'
            %Create and load the table on the current figure
            coordnew=[num2cell(true(length(LAT),1)) num2cell(LAT) num2cell(LON) num2cell(ELE) cellstr(utm_zone)];%Data
            cnames={'Selection','Lat(°)', 'Lon(°)','Elev(m)','Zones'};%Column names
            rnames=station';%Rownames
             dataCoord=uitable('Parent',Fig1,'Data',coordnew,'ColumnName',cnames,... 
            'Columnformat',columnformat,'RowName',rnames,'units','normalized','Position',[.37 .015 .285 .2]);
        case 'm'
            %Create and load the table on the current figure
            coordnew=[num2cell(true(length(xs),1)) num2cell(xs) num2cell(ys) num2cell(zs) cellstr(utm_zone)];%Data
            cnames={'Selection','X(km)', 'Y(km)','Z(km)','Zones'};%Column names
            rnames=station';%Rownames
            dataCoord=uitable('Parent',Fig1,'Data',coordnew,'ColumnName',cnames,... 
            'Columnformat',columnformat,'RowName',rnames,'units','normalized','Position',[.37 .015 .285 .2]);      
    end
    %Set the table as editable
    set(dataCoord,'ColumnEditable',true(1,4))
catch
    %Error assessment
    if isempty(params.coord) | isempty(station);set(text1, 'String','Error. No coordinate file directory or No correct coordinate file.');drawnow;return;end  
end

%% Other functions (nested)
%%%%Upload the new stations coordinates%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function upCoord(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    try
    coord_sys=params.coord_sys;
    %Upload the table on the basis of the geographic system (wgs84/Utm )
    switch coord_sys
        case 'degrees'
            data=get(dataCoord,'Data');%Get the new coordinates of the stations
            station=get(dataCoord,'RowName');%Get the stations names
            %Upload the data
            LAT=cell2mat(data(:,2));LON=cell2mat(data(:,3));ELE=cell2mat(data(:,4));idx=cell2mat(data(:,1));utm_zone=cell2mat(data(:,5));
            LAT=LAT(idx);LON=LON(idx);ELE=ELE(idx);station=station(idx);utm_zone=utm_zone(idx,:);
            %Command message
            set(text1, 'String', 'Station coordinates are uploaded.');
        case 'm'
            data=get(dataCoord,'Data');%Get the new coordinates of the stations
            station=get(dataCoord,'RowName');%Get the stations names
             %Upload the data
            xs=cell2mat(data(:,2));ys=cell2mat(data(:,3));zs=cell2mat(data(:,4));idx=cell2mat(data(:,1));utm_zone=cell2mat(data(:,5));
            xs=xs(idx);ys=ys(idx);zs=zs(idx);station=station(idx);utm_zone=utm_zone(idx,:);
            %Command message
            set(text1, 'String', 'Station coordinates are uploaded.');
    end
    catch
        %Error assessment
        if isempty(data);set(text1, 'String','Error. No coordinate file directory or No correct coordinate file.');drawnow;return;end  
    end
end

%==========================================================================
%%%%Save the new stations coordinates%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savCoord(~,~)
   text1.String=[];set(text1);drawnow;%Clear the Window message
   coord_sys=params.coord_sys;
   %Save the table on the basis of the geographic system (wgs84/Utm )
    switch coord_sys
        case 'm'
            %Error control
            if isempty(xs);set(text1, 'String','Error. No correct coordinate file.');drawnow;return;end
            %Sort the station coordinates
            [station,idx]=sort(station);
            xs=xs(idx);
            ys=ys(idx);
            zs=zs(idx);
            utm_zone=utm_zone(idx,:);
            %Save the station coordinates
            save(params.coord,'station','xs','ys','zs','utm_zone');
            %Convert WGS84 coordinates (Latitude, Longitude) into UTM coordinates
            [LAT,LON]=utm2deg(xs*1000,ys*1000,utm_zone);
            ELE=zs*1000;%Expressed in m 
            %Save the station coordinates
            save([params.coord(1:end-5) 'deg.mat'],'station','LAT','LON','ELE','utm_zone'); 
        case 'degrees'
            %Error control
            if isempty(LAT);set(text1, 'String','Error. No correct coordinate file.');drawnow;return;end
            %Sort the station coordinates
            [station,idx]=sort(station);
            LAT=LAT(idx);
            LON=LON(idx);
            ELE=ELE(idx);
            utm_zone=utm_zone(idx,:);
            %Save the station coordinates
            save(params.coord,'station','LAT','LON','ELE','utm_zone'); 
            %Convert WGS84 coordinates (Latitude, Longitude) into UTM coordinates
            [xs,ys,~]=wgs2utm(LAT,LON);
            zs=ELE;
            %Express coordinates in km 
            xs=xs./1000;
            ys=ys./1000;
            zs=zs./1000;
            %Save the station coordinates
            save([params.coord(1:end-7) 'm.mat'],'station','xs','ys','zs','utm_zone'); 
    end
    %Command window
    set(text1, 'String', 'Station coordinates are saved.');
end

%==========================================================================
%%%%%%%%%%%%%%%%Beam pattern analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calc(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes2)%Clear the Beam Pattern axes
    %Acquisition and validation of the analysis parameters
    fmin=str2num(get(h.fmin,'string'));if fmin<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    fmax=str2num(get(h.fmax,'string'));if fmax<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    step=str2num(get(h.step,'string'));if step<=0;set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    smin=str2num(get(h.smin,'string'));if smin==0;set(text1, 'String','Error. Invalid analysis slowness.');drawnow;return;end
    smax=str2num(get(h.smax,'string'));if smax==0;set(text1, 'String','Error. Invalid analysis slowness.');drawnow;return;end
    ns=str2num(get(h.ns,'string'));if ns<=0;set(text1, 'String','Error. Invalid analysis slowness.');drawnow;return;end
    F=fmin:step:fmax;if isempty(F);set(text1, 'String','Error. Invalid analysis frequencies.');drawnow;return;end
    %Command message
    textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
    set(text1, 'String', textLabel); 
    drawnow
    BeamP={};%Initialize the Beam Pattern matrix
    Sx=[];Sy=[];%Initialize the horizontal slowness vectors
    slidval=[];%Initialize the value of dynamic plot
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
    
    %Loop through the frequencies analysis
    for f=fmin:step:fmax
        %Apply the Beam Pattern algorithm
        [B,sx,sy]=BeamPatternS(xs',ys',smin,smax,f,ns);
        %Concatenate the results
        BeamP=[BeamP;B];
        Sx=[Sx sx'];
        Sy=[Sy sy'];
    end
    %Set the block of results
    slidval=length(F);
    %Plot the results in the beam pattern axes
    ploting2(slidval,axes2)
    %Command window
    set(text1, 'String','Calculation finished.');
    drawnow
end

%==========================================================================
%%%%%%%%Map construction and viewing the station coordinates%%%%%%%%%%%%%%%
function calc2(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    cla(axes1)%Clear the Map axes
    try
        %Acquisition and validation of the parameters
        zmin=str2num(get(h.zmin,'string'));if zmin<=0;set(text1, 'String','Error. Invalid elevation range.');drawnow;return;end
        zmax=str2num(get(h.zmax,'string'));if zmax<=0;set(text1, 'String','Error. Invalid elevation range.');drawnow;return;end
        zstep=str2num(get(h.zstep,'string'));if zstep<=0;set(text1, 'String','Error. Invalid elevation range.');drawnow;return;end
        if zmax<=zmin;set(text1, 'String','Error. Invalid elevation range.');drawnow;return;end  
        if isempty(params.coord) | isempty(station);set(text1, 'String','Error. No coordinate file directory or No correct coordinate file.');drawnow;return;end
        %Command message
        textLabel = sprintf('Calculation started. Be patient ! \n After finshing the result will be appear on the screen.');
        set(text1, 'String', textLabel); 
        drawnow
        %DEM construction and plot the station coordinates
        if ~isempty(params.map) 
            %Execute one of two groups of routines (coordinates in wgs84/UTM)
            switch params.coord_sys
                case 'm'
                    %Error control
                    if isempty(xs);set(text1, 'String','Error. No correct coordinate file.');drawnow;return;end
                    cla(axes1)%Clear the Map axes
                    load('mapfile_m.mat','X','Y','A')%Load the DEM variables (Utm system)
                    %Plot the DEM 
                    contourf(axes1,X,Y,A,(zmin:zstep:zmax),'k'); 
                    colorbar(axes1)
                    hold(axes1,'on')
                    view(axes1,[0 90])
                    %Plot the station coordinates and names (Utm system)
                    plot3(axes1,xs,ys,zs,'ok','MarkerFacecolor','r')
                    text(axes1,xs,ys,station,'VerticalAlignment','bottom','HorizontalAlignment','right')
                    hold(axes1,'off')
                    xlabel(axes1,'X (km)')
                    ylabel(axes1,'Y (km)')
                    set(axes1,'XLim',[xLimits(1)/1000 xLimits(2)/1000],'YLim',[yLimits(1)/1000 yLimits(2)/1000],'fontsize',8);
                    legend(axes1,{strcat(num2str(zstep*1000),'meters');'Stations'})
                    axis(axes1,'equal');
                    axis(axes1,'tight');
                case 'degrees'
                    %Error control
                    if isempty(LAT);set(text1, 'String','Error. No correct coordinate file.');drawnow;return;end
                    cla(axes1)%Clear the Map axes
                    load('mapfile_deg.mat','X','Y','A')%Load the DEM variables (wgs84 system)
                    %Plot the DEM 
                    contourf(axes1,X,Y,A*1000,(zmin*1000:zstep*1000:zmax*1000),'k')
                    colorbar(axes1)
                    hold(axes1,'on')
                    view(axes1,[0 90])
                    %Plot the station coordinates and names (wgs84 system)
                    plot3(axes1,LON,LAT,ELE,'ok','MarkerFacecolor','r')
                    text(axes1,LON,LAT,station,'VerticalAlignment','bottom','HorizontalAlignment','right')
                    hold(axes1,'off')
                    xlabel(axes1,'Longitude (degrees)')
                    ylabel(axes1,'Latitude (degrees)')
                    set(axes1,'XLim',[xLimits(1) xLimits(2)],'YLim',[yLimits(1) yLimits(2)],'fontsize',8);
                    legend(axes1,{strcat(num2str(zstep*1000),'meters');'Stations'})
                    axis(axes1,'equal');
                    axis(axes1,'tight');
            end
        %Plot only the station coordinates    
        elseif isempty(params.map)
            %Execute one of two groups of routines (coordinates in wgs84/UTM)
            switch params.coord_sys
                case 'degrees'
                    %Error control
                    if isempty(LAT);set(text1, 'String','Error. No correct coordinate file.');drawnow;return;end
                    cla(axes1)%Clear the Map axes
                    %Plot the station coordinates and names (wgs84 system)
                    plot(axes1,LON,LAT,'ko','MarkerFaceColor','r')
                    hold(axes1,'on')
                    text(axes1,LON,LAT,station,'VerticalAlignment','bottom','HorizontalAlignment','right')
                    hold(axes1,'off')
                    grid(axes1,'on')
                    grid(axes1,'minor')
                    xlabel(axes1,'Latitude (decimal degrees)')
                    ylabel(axes1,'Longitude (decimal degrees)')
                    set(axes1,'XLim',[min(LON)-0.001 max(LON)+0.001],'YLim',[min(LAT)-0.001 max(LAT)+0.001],'fontsize',8);
                    legend(axes1,['Stations'])
                    axis(axes1,'equal');
                    axis(axes1,'tight');
                case 'm'
                    %Error control
                    if isempty(xs);set(text1, 'String','Error. No correct coordinate file.');drawnow;return;end
                    cla(axes1)%Clear the Map axes
                    %Plot the station coordinates and names (Utm system)
                    plot(axes1,xs,ys,'ko','MarkerFaceColor','r') 
                    hold(axes1,'on')
                    text(axes1,xs,ys,station,'VerticalAlignment','bottom','HorizontalAlignment','right')
                    hold(axes1,'off')
                    grid(axes1,'on')
                    grid(axes1,'minor')
                    xlabel(axes1,'X (km)')
                    ylabel(axes1,'Y (km)')
                    set(axes1,'XLim',[min(xs)-0.1 max(xs)+0.1],'YLim',[min(ys)-0.1 max(ys)+0.1],'fontsize',8);
                    legend(axes1,['Stations'])
                    axis(axes1,'equal');
                    axis(axes1,'tight');
            end
        end
        %Command window
        set(text1, 'String','Calculation finished.');
        drawnow
    catch
        %Error assessment
        if isempty(X) | isempty(Y) | isempty(A);set(text1, 'String','Error. No map file directory or No correct map file.');drawnow;return;end
    end
end

%==========================================================================
%%%%%%%%Plot the Beam pattern results%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploting2(ii,ax)
    hold(ax,'on')
    cla(ax)%Clear the Beam Pattern axes   
    colormap(ax,parula)
    colorbar(ax)
    %Plot the results of the analysis 
    pcolor(ax,Sx(:,ii),Sy(:,ii),BeamP{ii});shading(ax,'interp');
    contour(ax,Sx(:,ii),Sy(:,ii),BeamP{ii}); 
    title(ax,['Freq=' num2str(F(ii)) ' Hz'])
    xlabel(ax,'S_x (s/km)', 'fontsize', 7);
    ylabel(ax,'S_y (s/km)', 'fontsize', 7);
    set(ax,'XLim',[smin smax],'YLim',[smin smax],'fontsize',7);
    axis(ax,'equal')
    axis(ax,'square')
    hold(ax,'off')
    
end

%==========================================================================
%%%%Refresh the Beam Pattern axes through the pre button%%%%%%%%%%%%%%%%%%%
function beam1(~,~)
     %Error control
     if isempty(cell2mat(BeamP));set(text1, 'String','Error. Beam Pattern matrix is empty.');drawnow;return;end
     %Decrease the slidval value
     slidval=slidval-1;if slidval<=0;slidval=1;end
     %Plot the new block of Beam Pattern parameters
     ploting2(slidval,axes2);        
end

%==========================================================================
%%%%Refresh the Beam Pattern axes through the next button%%%%%%%%%%%%%%%%%%
function beam2(~,~)
    %Error control
    if isempty(cell2mat(BeamP));set(text1, 'String','Error. Beam Pattern matrix is empty.');drawnow;return;end
    %Increase the slidval value
    slidval=slidval+1;if slidval>length(F);slidval=length(F);end
    %Plot the new block of Beam Pattern parameters
    ploting2(slidval,axes2);        
end

%==========================================================================
%%%%Control and reset automatically analysis parameters%%%%%%%%%%%%%%%%%%%%
function deflimits(~,~)
    %Acquisition analysis parameters
    fmin=str2num(get(h.fmin,'string'));
    fmax=str2num(get(h.fmax,'string'));
    step=str2num(get(h.step,'string'));
    smin=str2num(get(h.smin,'string'));
    smax=str2num(get(h.smax,'string'));
    ns=str2num(get(h.ns,'string'));
    zmin=str2num(get(h.zmin,'string'));
    zmax=str2num(get(h.zmax,'string'));
    zstep=str2num(get(h.zstep,'string'));
    
    %Control and eventually reset parameters
    if (size(fmin,1)==0);set(h.fmin,'string',params.fmin);end
    if (size(fmax,1)==0);set(h.fmax,'string',params.fmax);end
    if (size(step,1)==0);set(h.step,'string',params.step);end
    if (size(smin,1)==0);set(h.smin,'string',params.smin);end
    if (size(smax,1)==0);set(h.smax,'string',params.smax);end
    if (size(ns,1)==0);set(h.ns,'string',params.ns);end
    if (size(zmin,1)==0);set(h.zmin,'string','1.0');end
    if (size(zmax,1)==0);set(h.zmax,'string','3.3');end
    if (size(zstep,1)==0);set(h.zstep,'string','0.05');end
end

%==========================================================================
%%%%Save the Beam Pattern plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save1(~,~)
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %Error control
    if isempty(cell2mat(BeamP));set(text1, 'String','Error. Beam Pattern matrix is empty.');drawnow;return;end
    ii=slidval;
    %Plot the results of the analysis 
    new_figure=figure('Name','Beam Pattern Plot','NumberTitle','off');   
    colormap(parula)
    pcolor(Sx(:,ii),Sy(:,ii),BeamP{ii});shading interp; hold on
    contour(Sx(:,ii),Sy(:,ii),BeamP{ii}); 
    title(['Freq=' num2str(F(ii)) ' Hz'])
    xlabel('S_x (s/km)', 'fontsize', 7);
    ylabel('S_y (s/km)', 'fontsize', 7);
    set(gca,'XLim',[smin smax],'YLim',[smin smax],'fontsize',12);
    axis equal
    axis square
    colorbar
end

%==========================================================================
%%%%Save the Map plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save2(~,~)   
    text1.String=[];set(text1);drawnow;%Clear the Window message
    %DEM construction and plot the station coordinates
    if ~isempty(params.map)
        %Execute one of two groups of routines (coordinates in wgs84/UTM)
        switch params.coord_sys
            case 'm'
                %Error control
                if isempty(X) | isempty(xs);set(text1, 'String','Error. Map matrix is empty.');drawnow;return;end
                new_figure=figure('Name','Station Map','NumberTitle','off');
                %Plot the DEM (Utm system) 
                contourf(X,Y,A,(zmin:zstep:zmax),'k'); 
                colorbar
                hold('on')
                view([0 90])
                %Plot the station coordinates and names (Utm system)
                plot3(xs,ys,zs,'ok','MarkerFacecolor','r')
                text(xs,ys,station,'VerticalAlignment','bottom','HorizontalAlignment','right')
                hold('off')
                xlabel('X (km)')
                ylabel('Y (km)')
                set(gca,'XLim',[xLimits(1)/1000 xLimits(2)/1000],'YLim',[yLimits(1)/1000 yLimits(2)/1000],'fontsize',8);
                legend(['50meters';'Stations'])
                axis('equal');
                axis('tight');
            case 'degrees'
                %Error control
                if isempty(X) | isempty(LAT);set(text1, 'String','Error. Map matrix is empty.');drawnow;return;end
                new_figure=figure('Name','Station Map','NumberTitle','off');
                %Plot the DEM (wgs84 system)  
                contourf(X,Y,A*1000,(zmin*1000:zstep*1000:zmax*1000),'k'); 
                colorbar
                hold('on')
                view([0 90])
                %Plot the station coordinates and names (wgs84 system)
                plot3(LON,LAT,ELE,'ok','MarkerFacecolor','r')
                text(LON,LAT,station,'VerticalAlignment','bottom','HorizontalAlignment','right')
                hold('off')
                xlabel('Longitude (degrees)')
                ylabel('Latitude (degrees)')
                set(gca,'XLim',[xLimits(1) xLimits(2)],'YLim',[yLimits(1) yLimits(2)],'fontsize',8);
                legend(['50meters';'Stations'])  
                axis('equal');
                axis('tight');
        end
    %Plot only the station coordinates  
    elseif isempty(params.map)
        %Execute one of two groups of routines (coordinates in wgs84/UTM)
        switch params.coord_sys
            case 'degrees'
                %Error control
                if  isempty(LAT);set(text1, 'String','Error. Map matrix is empty.');drawnow;return;end
                new_figure=figure('Name','Station Map','NumberTitle','off');
                %Plot the station coordinates and names (wgs84 system)
                plot(LON,LAT,'ko','MarkerFaceColor','r')
                hold('on')
                text(LON,LAT,station,'VerticalAlignment','bottom','HorizontalAlignment','right')
                hold('off')
                grid('on')
                grid('minor')
                xlabel('Latitude (decimal degrees)')
                ylabel('Longitude (decimal degrees)')
                set(gca,'XLim',[min(LON)-0.001 max(LON)+0.001],'YLim',[min(LAT)-0.001 max(LAT)+0.001],'fontsize',8);
                legend('Stations')
                axis('equal');
                axis('tight');
            case'm'
                %Error control
                if isempty(X) & isempty(xs);set(text1, 'String','Error. Map matrix is empty.');drawnow;return;end
                new_figure=figure('Name','Station Map','NumberTitle','off');
                %Plot the station coordinates and names (Utm system)
                plot(xs,ys,'ko','MarkerFaceColor','r') 
                hold('on')
                text(xs,ys,station,'VerticalAlignment','bottom','HorizontalAlignment','right')
                hold('off')
                grid('on')
                grid('minor')
                xlabel('X (km)')
                ylabel('Y (km)')
                set(gca,'XLim',[min(xs)-0.1 max(xs)+0.1],'YLim',[min(ys)-0.1 max(ys)+0.1],'fontsize',8);
                legend('Stations')
                axis('equal');
                axis('tight');
        end
    end
end
end