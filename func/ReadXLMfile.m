%%%%%%%%%%%%% Read XML file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% November 2022 
% Vittorio Minio
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xlm: XML file directory+name
% T: Matlab struct related to XML file
% lat, lon, ele: coordinates in wgs84 system
% poles: poles of the transfer function
% zers, zeros: zeros of the transfer function
% k: calibration coefficient (rad/s)
% vvel: coefficient for conversion from Volts to m/s (V s/m)
% bmw: coefficient for conversion from counts to Volts (V/counts)
%az,dip,depth: sensor orientation 

function [sta,ch,lat,lon,ele,vvel,bmw,k,az,dip,depth,units,poles,zers]=ReadXLMfile(xlm)
    %Read Matlab struct
    T=readstruct(xlm);
    %Extract all information  
    sta=T.Network.Station.codeAttribute;
    ch=T.Network.Station.Channel.codeAttribute;
    lat=T.Network.Station.Channel.Latitude.Text;
    lon=T.Network.Station.Channel.Longitude.Text;
    ele=T.Network.Station.Channel.Elevation.Text;
    sensivity=T.Network.Station.Channel.Response.InstrumentSensitivity.Value;
    vvel=T.Network.Station.Channel.Response.Stage(1).StageGain.Value;
    bmw= vvel/sensivity;
    k=T.Network.Station.Channel.Response.Stage(1).PolesZeros.NormalizationFactor;
    zer_struct=T.Network.Station.Channel.Response.Stage(1).PolesZeros.Zero;nzer=length(zer_struct);
    pol_struct=T.Network.Station.Channel.Response.Stage(1).PolesZeros.Pole;npol=length(pol_struct);
    poles=complex(zeros(npol,1),0);
    dip=T.Network.Station.Channel.Dip.Text;
    az=T.Network.Station.Channel.Azimuth.Text;
    depth=T.Network.Station.Channel.Depth.Text;
    units=T.Network.Station.Channel.Response.InstrumentSensitivity.InputUnits.Name;
    for pl=1:npol
        poles(pl)=complex(pol_struct(pl).Real.Text,pol_struct(pl).Imaginary.Text);
    end
    zers=complex(zeros(nzer,1),0);
    for zr=1:nzer
        zers(zr)=complex(zer_struct(zr).Real.Text,zer_struct(zr).Imaginary.Text);
    end
end