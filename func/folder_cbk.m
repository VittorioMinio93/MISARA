%%%%%%%%%%%%% Set user paths and folders %%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h2: UI control parameter
% X,Y,A: georeferenced DEM grid NxM (UTM system/WGS84 system)
% R: spatial referencing object
% xLimits, yLimits: georeferenced limits of the DEM grid in m/decimal degrees
% utm_zone:  UTM zone of reference 
% dem_UTM, ydem_UTM: coordinate vectors of the DEM grid (km) 
% xx, yy: coordinate vectors of the DEM grid (decimal degrees)
% station: name of the stations vector
% xs, ys, zs: coordinate of the stations (UTM system)
% LAT, LON, ELE: coordinate of the stations (WGS84 system)
% Data: table of the instrumental response correction parameters 
%%%%%%%%%%%h2.variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data: the Main Input folder
% save: the Main Output folder
% sacin: the Input folder for the pre-existing repository (Create dataset module)
% savein: the Output folder for the pre-existing repository (Create dataset module)
% irisout: the Output folder for the Iris repository
% map: the path of the DEM file (.tif)
% coord: the path of the coordinate file (.mat) 
% poles: the path of the instrumental response correction parameters (.mat)
% data2: the Input folder for the Signal viewer module
% save2: the Output folder for the Signal viewer module
function folder_cbk(hobj,~)
%% Declaration global parameters
global h2 comp text1
%%Execute one of ten groups of routines
%% Set user path and folders
switch get(hobj,'Tag')
%--- Main Input
case 'data_folder'      
    tmp= uigetdir(get(h2.data,'String'));
    if tmp~=0
        set(h2.data,'String',tmp)
    end    
    
    
%--- Main Output
case 'save_folder'      
    tmp= uigetdir(get(h2.save,'String'));
    if tmp~=0
        set(h2.save,'String',tmp)
    end 
    
%--- Pre-existing files Output
case 'savein'      
    tmp= uigetdir(get(h2.savein,'String'));
    if tmp~=0
        set(h2.savein,'String',tmp)
    end  
    
%--- Pre-existing files Input
case 'sacin'      
    tmp= uigetdir(get(h2.sacin,'String'));
    if tmp~=0
        set(h2.sacin,'String',tmp)
    end 
%--- Pre-existing files Input
case 'pathsta'      
    tmp= uigetdir(get(h2.pathsta,'String'));
    if tmp~=0
        set(h2.pathsta,'String',tmp)
    end
%--- .XLM files Input
case 'xlm'      
    [filename, pathname] = uigetfile( ...
        {[h2.staz_reff.String '.*' comp '*.xml'],'Data Files (*.xml)';...
        '*.*',  'All Files (*.*)'}, 'Pick a file',get(h2.xlm,'String'));%Select the XMLfile
    if filename~=0
        xlm=[pathname,filename];
        set(h2.xlm,'String',xlm)
        %%Read XML file
        [sta,ch,lat,lon,ele,vvel,bmw,k,az,dip,depth,units,poles,zers]=ReadXLMfile(xlm);       
        
        set(h2.lat_coord,'String',lat)
        set(h2.lon_coord,'String',lon)
        set(h2.ele_coord,'String',ele)
        set(h2.vvel,'String',vvel)
        set(h2.bmw,'String',bmw)
        set(h2.k,'String', k)
        set(h2.pols,'String',num2str(poles))
        set(h2.zers,'String',num2str(zers))

        h2.dip=dip;
        h2.az=az;
        h2.depth=depth;
        h2.units=units;
        %Command message
        textLabel = sprintf('XLM file is read.');
        set(text1, 'String', textLabel);
    else
        %Error message
        textLabel = sprintf('Error. No XLM files are read.');
        set(text1, 'String', textLabel);
    end 
    
%--- Iris files Output
case 'irisout'      
    tmp= uigetdir(get(h2.irisout,'String'));
    if tmp~=0
        set(h2.irisout,'String',tmp)
    end
    
%--- Signal viewer Input
case 'data2_folder'      
    tmp= uigetdir(get(h2.data2,'String'));
    if tmp~=0
        set(h2.data2,'String',tmp)
    end 
    
%--- Signal viewer Output
case 'save2_folder'      
    tmp= uigetdir(get(h2.save2,'String'));
    if tmp~=0
        set(h2.save2,'String',tmp)
    end 
%% Set user path and folders/ Sort and save analysis parameters
%--- DEM file
case 'map'
    [filename, pathname] = uigetfile( ...
        {'*.tif','Data Files (*.tif)';...
        '*.*',  'All Files (*.*)'}, 'Pick a file',get(h2.map,'String'));%Select the DEM file
    if filename~=0
        set(h2.map,'String',[pathname,filename])%Set the path of the DEM file
        map_file=get(h2.map,'string');
        [A, R] = geotiffread(map_file);%Read the DEM file selected
        A = double(A)/1000;%Expressed in km                 
       [N,M]= size(A); %Get the dimensions of the elevetion grid NxM
       xLimits=R.LongitudeLimits;%Expressed in decimal degrees 
       yLimits=R.LatitudeLimits;%Expressed in decimal degrees 
       yy= linspace(yLimits(2), yLimits(1), N);      
       xx= linspace(xLimits(1), xLimits(2), M);
      [X,Y]=meshgrid(xx,yy);%Get georeferenced DEM grid NxM (decimal degrees)
      [~,~,utm_zone] = deg2utm(R.LatitudeLimits,R.LongitudeLimits);%Return the UTM zone fo reference
       save(strcat(cd,'/parameters/mapfile_deg.mat'),'X','Y','A','xLimits','yLimits','utm_zone')%Save the DEM file variables (WGS84)
       [xLimits,yLimits,utm_zone] = deg2utm(R.LatitudeLimits,R.LongitudeLimits);%Convert vectors of Lat/Lon (WGS84) coordinates into UTM vectors
       ydem_UTM= linspace(yLimits(2), yLimits(1), N)/1000;
       xdem_UTM= linspace(xLimits(1), xLimits(2), M)/1000;
       [X,Y]=meshgrid(xdem_UTM,ydem_UTM);%Get georeferenced DEM grid NxM (km)
       save(strcat(cd,'/parameters/mapfile_m.mat'),'X','Y','A','xLimits','yLimits','utm_zone')%Save the DEM file variables (UTM system)
    end
    
%--- Coordinate file
case 'coord'
    [filename, pathname] = uigetfile( ...
        {'*.mat;*.dat','Data Files (*.mat)';...
        '*.*',  'All Files (*.*)'}, 'Pick a file',get(h2.coord,'String'));%Select the coordinate file of the stations
    if filename~=0
        set(h2.coord,'String',[pathname,filename])%Set the path of the coordinate file
        coord=get(h2.coord,'string');
        xs=[];ys=[];zs=[];station=[];LAT=[];LON=[];ELE=[];%Initialize the coordinate station vectors
        load(coord)%Load the coordinate file selected
        %Execute one of two groups of routines based on the Geographic system (UTM/WGS84)
        if  ~isempty(xs)
            %Sort coordinate station vector (UTM system)
            [station,idx]=sort(station);
            xs=xs(idx);
            ys=ys(idx);
            zs=zs(idx);
            utm_zone=utm_zone(idx,:);
            %Save and clear the coordinate station vectors (UTM system) 
            save(coord,'station','xs','ys','zs','utm_zone');
            clear station xs ys zs
        elseif ~isempty(LAT)
            %Sort coordinate station vector (WGS84 system)
            [station,idx]=sort(station);
            LAT=LAT(idx);
            LON=LON(idx);
            ELE=ELE(idx);
            utm_zone=utm_zone(idx,:);
            %Save and clear the coordinate station vectors (WGS84 system)
            save(coord,'station','LAT','LON','ELE','utm_zone');  
            clear station LAT LON ELE
        end
    end
    
%--- Instrumental response correction parameters
case 'poles'
    [filename, pathname] = uigetfile( ...
        {'*.mat','Data Files (*.mat)';...
        '*.*',  'All Files (*.*)'}, 'Pick a file',get(h2.poles,'String'));%Select the instrumental response correction parameters file
    if filename~=0
        set(h2.poles,'String',[pathname,filename])%Set the path of the instrumental response correction parameters file
        poles=get(h2.poles,'string');
        load(poles)%Load the instrumental response correction parameters file 
        Data=sortrows(Data,1,'ascend');%Sort the table of the instrumental response correction parameters
        %Save and clear the table of the instrumental response correction parameters
        save(poles,'Data');
        clear Data
    end 
end
end