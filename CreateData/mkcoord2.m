%%%%%%%%%%%%% Save station coordinate strucuters%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pathout: the Output folder 
% pathcoord: the Output coordinates directory
% sta: the name of the reference station 
% latitude/LAT, longitude/LON, elevation/ELE: coordinates in wgs84 system
% xs, ys, zs, utm_zone: coordinates in Utm system 

function mkcoord2(pathout,sta,latitude,longitude,elevation)
%Create the Output directory
pathcoord=strcat(pathout,'/coordinate/');if ~exist(pathcoord, 'dir');mkdir(pathcoord);end
%Search the station coordinates file
p=dir(strcat(pathcoord,'coordinate_deg.mat'));
%Control and eventually set the station coordinates
if isempty(latitude);latitude=NaN;end
if isempty(longitude);longitude=NaN;end
if isempty(elevation);elevation=NaN;end
if length(p)==0
    %Create the new station coordinates vectors 
    LAT=latitude;LON=longitude;ELE=elevation;station={sta};
else
    %Load and Upload the new station coordinates values 
    load(strcat(pathcoord,'coordinate_deg.mat'))
    staidx=strcmp(sta,station);
    LAT(staidx)=[];LON(staidx)=[];ELE(staidx)=[];station(staidx)=[];
    LAT=[LAT;latitude];LON=[LON;longitude];ELE=[ELE;elevation];station=[station;{sta}];
end
%Delete the no unique values
[~,idx]=unique(station);
LAT=LAT(idx);LON=LON(idx);ELE=ELE(idx);station=station(idx);
%Define the Utm zone
if isnan(LAT) | isnan(LON) 
     utm_zone=NaN;
 elseif ~isnan(LAT) & ~isnan(LON) 
 [~,~,utm_zone] = deg2utm(LAT,LON);
end
%Save the station coordinates vectors (wgs84 system)
save([pathcoord 'coordinate_deg.mat'],'LAT','LON','ELE','station','utm_zone')
%Convert WGS84 coordinates (Latitude, Longitude) into UTM coordinates
 if isnan(LAT) | isnan(LON) 
     xs=LAT;ys=LON;
 elseif ~isnan(LAT) & ~isnan(LON) 
 [xs,ys,utm_zone] = deg2utm(LAT,LON);
 end
%Express coordinates in km 
xs=xs./1000;ys=ys./1000;zs=ELE./1000;
%Save the station coordinates vectors (Utm system) 
save([pathcoord 'coordinate_m.mat'],'xs','ys','zs','station','utm_zone')
end