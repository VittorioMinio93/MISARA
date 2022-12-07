%%%%%%%%%%%%% Get and save user parameters in the Home window%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%params: user parameters
%h2: UI control parameters
%%%%%%%%%%%h2.variables/params.variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% radiobutton: radio button within a button group 
% radio: button group 
%-------------------------------- General settings
%--- General info
% coord: the path of the coordinate file (.mat) 
% map: the path of the DEM file (.tif)
% sta_ref: the name of the reference station 
% component: the name of the component/channel
% comp_chan: the station system used (component/channel)
%--- Data file
% data: the Input folder
% save: the Output folder
% outputformat: the format (.mat/.txt) of the output files 
%-------------------------------- Instrumental settings
%--- Coordinate system/Beam Pattern setting
% coord_sys: the coordinate system (metrics/degrees)
% fmin: the minimum frequency for the Beam pattern analysis (Hz)
% fmax: the maximum frequency for the Beam pattern analysis(Hz)
% step: the frequency step for the Beam pattern analysis (Hz)
% smin: the minimum slowness for the Beam pattern analysis (s/km)
% smax: the miximum slowness for the Beam pattern analysis (s/km)
% ns: the dimensions of the square grid slowness for the Beam pattern analysis
%--- Instrument correction setting
% poles: the path of the instrumental response correction parameters (.mat)
% data2: the Input folder
% save2: the Output folder
%-------------------------------- Preliminary analysis settings
%--- Spectral setting
% w_spec: the analysis window for the Spectrogram (s)
% np: the dimension of the spectrums
% f_high_pass: the frequency for the high-pass filter (Hz)
% n_sta: the maximum number of stations to use for the Spectral coherence
% spec_hour: the number of spectrums for each trace
%--- RMS setting
% w_rms: the analysis window for the RMS analysis (s)
% f1: the minimum frequency for the RMS analysis (Hz)
% f2: the maximum frequency for the RMS analysis (Hz)
% step_rms: the frequency step for the RMS analysis (Hz)
% factor: the temporal factor for the moving average of the RMS (hours) 
%--- Polarization setting
% w_pol: the analysis window for the polarization analysis (s)
% f1_pol: the minimum frequency for the polarization analysis (Hz)
% f2_pol: the maximum frequency for the polarization analysis (Hz)
% step_pol: the frequency step for the polarization analysis (Hz)
% factor_pol: the temporal factor for the moving average of the polarization attributes (hours)
%--- Detection setting
% wshort: the short-term window or the pre-trigger window (s)
% wlong: the long-term window or the post-trigger window (s)
% f1_det: the minimum frequency for the band-pass filter (Hz)
% f2_det: the maximum frequency for the band-pass filter (Hz)
% treshold: the threshold value of the detection
% f1_down:  the minimum frequency for the band-pass filter (Only Salped) (Hz)
% f2_down: the maximum frequency for the band-pass filter (Only Salped) (Hz)
% f1_up: the minimum frequency for the band-pass filter (Only Salped) (Hz)
% f2_up: the maximum frequency for the band-pass filter (Only Salped) (Hz)
%-------------------------------- Array analysis settings
%--- ZLC setting
% w_zlc: the analysis window for the ZLC analysis (s)
% f1_zlc: the minimum frequency for the ZLC analysis (Hz)
% f2_zlc: the maximum frequency for the ZLC analysis (Hz)
% step_zlc: the frequency step for the ZLC analysis (Hz)
% vel_zlc: the velocity value (km/s)
% splint: the Spline interpolation value (1 o 0)
% maxlags: the max lags for the Correlation coefficient calculation (s)
%--- Radial Semblance/Semblance setting
% w_sem: the analysis window for the Radial Semblance/Semblance analysis (s)
% f0: the central frequency for the seismic amplitude decays with distance (Hz)
% Q: the quality factor for the seismic amplitude decays with distance 
% vel_sem: the velocity value (km/s)
% n_sem: the exponent value 
% step_sem: the grid step size (km)
% topo_corr: the topographic correction (yes o no)
% freq_sem: the  frequency range for the band-pass filter (Hz)
% xmin: the minimum X-limit of the grid (km)
% xmax: the maximum X-limit of the grid (km)
% ymin: the minimum Y-limit of the grid (km)
% ymax: the maximum Y-limit of the grid (km)
% zmin: the minimum Z-limit of the grid (km)
% zmax: the maximum Z-limit of the grid (km)
%--- MUSIC setting
% freq0: the central frequency of the band of interest (Hz)
% nfreq: the number of frequencies of the band of interest
% m: the hamming value
% R: the R value
% pmax: the limit of the grid slowness (s/km) 
% pinc: the size step of the grid slowness (s/km) 
% ninicio: the time origin of the analysis window (s)  
% lwin: the analysis window for the MUSIC analysis (s) 
% advance: the advancement of the analysis window for the MUSIC analysis (%) 
% nsig: the number of source to use exclusively for the plotting 

function save_cbk(hobj,~)
%% Declaration global parameters
global params h2
%% Save user parameters
%%Execute one of four groups of settings
switch get(hobj,'Tag')
case 'save1'
%-------------------------------- General settings
%--- General info
    params.coord=deblank(get(h2.coord,'string'));
    params.map= deblank(get(h2.map,'string'));
    params.sta_ref= deblank(get(h2.sta_ref,'string'));
    
    tmp= get(h2.radiobutton(3),'selectedobject');%%Get the name of the component/channel from radiobutton     
    params.component= tmp.String; 
    
    tmp= get(h2.radiobutton(4),'selectedobject');%%Get the station system used (component/channel) from radiobutton     
    params.comp_chan= tmp.String;
    
%--- Data file 
    params.data= deblank(get(h2.data,'string'));
    params.save= deblank(get(h2.save,'string'));
    
    tmp= get(h2.radiobutton(1),'selectedobject');%%Get the format (.mat/.txt) of the output files from radiobutton    
    params.outputformat= tmp.String;
    
    %%Save User values
    save('./parameters/user_values.mat','params');%%Save the user values for the General settings  
    
case 'save2' 
%-------------------------------- Instrumental settings
%--- Coordinate system/Beam Pattern setting
    tmp= get(h2.radiobutton(2),'selectedobject');%%Get the coordinate system (UTM/WGS84) from radiobutton    
    params.coord_sys= tmp.String;
    params.fmin= str2num( get(h2.fmin,'string'));
    params.fmax= str2num( get(h2.fmax,'string'));
    params.step= str2num( get(h2.step,'string'));
    params.smin= str2num( get(h2.smin,'string'));
    params.smax= str2num( get(h2.smax,'string'));
    params.ns= str2num( get(h2.ns,'string'));

%--- Instrument correction setting    
    params.poles=deblank(get(h2.poles,'string'));
    params.data2= deblank(get(h2.data2,'string'));
    params.save2= deblank(get(h2.save2,'string'));
    
    %%Save User values
    save('./parameters/user_values.mat','params');%%Save the user values for the Instrumental settings  
    
    
case 'save3' 
%-------------------------------- Preliminary analysis settings
%--- Spectral setting
    params.w_spec= str2num( get(h2.w_spec,'string'));
    params.np= str2num( get(h2.np,'string'));
    params.f_high_pass= str2num( get(h2.f_high_pass,'string'));
    params.n_sta= str2num( get(h2.n_sta,'string'));
    params.spec_hour= str2num( get(h2.spec_hour,'string'));
%--- RMS setting
    params.w_rms= str2num( get(h2.w_rms,'string'));
    params.f1= str2num( get(h2.f1,'string'));
    params.f2= str2num( get(h2.f2,'string'));
    params.step_rms= str2num( get(h2.step_rms,'string'));
    params.factor= str2num( get(h2.factor,'string'));
%--- Polarization setting
    params.w_pol= str2num( get(h2.w_pol,'string'));
    params.f1_pol= str2num( get(h2.f1_pol,'string'));
    params.f2_pol= str2num( get(h2.f2_pol,'string'));
    params.step_pol= str2num( get(h2.step_pol,'string'));
    params.factor_pol= str2num( get(h2.factor_pol,'string'));
%--- Detection setting
    params.wshort=str2num( get(h2.wshort,'string'));
    params.wlong=str2num( get(h2.wlong,'string'));
    params.f1_det= str2num( get(h2.f1_det,'string'));
    params.f2_det= str2num( get(h2.f2_det,'string'));
    params.treshold= str2num( get(h2.treshold,'string'));
    params.f1_down= str2num( get(h2.f1_down,'string'));
    params.f2_down= str2num( get(h2.f2_down,'string'));
    params.f1_up= str2num( get(h2.f1_up,'string'));
    params.f2_up= str2num( get(h2.f2_up,'string'));
    
    %%Save User values
    save('./parameters/user_values.mat','params');%%Save the user values for the Preliminary analysis settings  
    
case 'save4'
%-------------------------------- Array analysis settings
%--- ZLC setting
    params.w_zlc= str2num( get(h2.w_zlc,'string'));
    params.f1_zlc= str2num( get(h2.f1_zlc,'string'));
    params.f2_zlc= str2num( get(h2.f2_zlc,'string'));
    params.step_zlc= str2num( get(h2.step_zlc,'string'));
    params.vel_zlc= str2num( get(h2.vel_zlc,'string'));
    params.splint= str2num( get(h2.splint,'string'));
    params.maxlags= str2num( get(h2.maxlags,'string'));
%--- Radial Semblance/Semblance setting
    params.w_sem= str2num( get(h2.w_sem,'string'));
    params.f0= str2num( get(h2.f0,'string'));
    params.Q= str2num( get(h2.Q,'string'));
    params.vel_sem= str2num( get(h2.vel_sem,'string'));
    params.n_sem= str2num( get(h2.n_sem,'string'));
    params.step_sem= str2num( get(h2.step_sem,'string'));
    params.topo_corr=get(h2.topo_corr,'string');
    params.freq_sem= deblank(get(h2.freq_sem,'string'));
    params.xmin= str2num( get(h2.xmin,'string'));
    params.xmax= str2num( get(h2.xmax,'string'));
    params.ymin= str2num( get(h2.ymin,'string'));
    params.ymax= str2num( get(h2.ymax,'string'));
    params.zmin= str2num( get(h2.zmin,'string'));
    params.zmax= str2num( get(h2.zmax,'string'));
%--- MUSIC setting
    params.freq0=str2num( get(h2.freq0,'string'));
    params.nfreq=str2num( get(h2.nfreq,'string'));
    params.m= str2num( get(h2.m,'string'));
    params.R= str2num( get(h2.R,'string'));
    params.pmax= str2num( get(h2.pmax,'string'));
    params.pinc= str2num( get(h2.pinc,'string'));
    params.ninicio= str2num( get(h2.ninicio,'string')); 
    params.lwin= str2num( get(h2.lwin,'string')); 
    params.advance= str2num( get(h2.advance,'string'));   
    params.nsig= str2num( get(h2.nsig,'string'));    
    
    %%Save User values
    save('./parameters/user_values.mat','params');%%Save the user values for the Array analysis settings
end
end