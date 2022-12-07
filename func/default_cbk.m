%%%%%%%%%%%%% Load and set default parameters in the Home window%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio

%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tmp: user parameters
% h2: UI control parameters
%%%%%%%%%%%h2.variables/tmp.params.variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function default_cbk(hobj,~)
%% Declaration global parameters
global h2

%% load default parameters
tmp= load('default_values.mat','params');

%% Set default parameters
%%Execute one of four groups of settings
switch get(hobj,'Tag')
case 'default1' 
%-------------------------------- General settings
%--- General info
    set(h2.coord,'string',tmp.params.coord);
    set(h2.map,'string',tmp.params.map);%%Set 
    set(h2.sta_ref,'string',tmp.params.sta_ref); 
    
    if strcmp(tmp.params.component,'Z')==1
        set(h2.radiobutton(3),'selectedobject',h2.radio(5));%%Set Z as the component/channel of reference  
    elseif strcmp(tmp.params.component,'N')==1
        set(h2.radiobutton(3),'selectedobject',h2.radio(6));%%Set N as the component/channel of reference
    elseif strcmp(tmp.params.component,'E')==1
        set(h2.radiobutton(3),'selectedobject',h2.radio(7));%%Set E asthe component/channel of reference 
    elseif strcmp(tmp.params.component,'F')==1
        set(h2.radiobutton(3),'selectedobject',h2.radio(10));%%Set F as the component/channel of reference
    end
    
   if strcmp(tmp.params.comp_chan,'Comp')==1
        set(h2.radiobutton(4),'selectedobject',h2.radio(8));%%Set Component as the station system of reference 
    elseif strcmp(tmp.params.comp_chan,'Chan')==1
        set(h2.radiobutton(4),'selectedobject',h2.radio(9));%%Set Channel as the station system of reference    
   end 
    
%--- Data file 
    set(h2.data,'string',tmp.params.data);
    set(h2.save,'string',tmp.params.save);
        
    if strcmp(tmp.params.outputformat,'.mat')==1
        set(h2.radiobutton(1),'selectedobject',h2.radio(1));%%Set .mat as the output format of reference     
    elseif strcmp(tmp.params.outputformat,'txt')==1
        set(h2.radiobutton(1),'selectedobject',h2.radio(2));%%Set .txt as the output format of reference     
    end 
    
case 'default2'
%-------------------------------- Instrumental settings
%--- Coordinate system/Beam Pattern setting
    if strcmp(tmp.params.coord_sys,'m')==1
        set(h2.radiobutton(2),'selectedobject',h2.radio(3)); %%Set UTM as the geographic system of reference  
    elseif strcmp(tmp.params.coord_sys,'degrees')==1
        set(h2.radiobutton(2),'selectedobject',h2.radio(4));%%Set WGS84 as the geographic system of reference 
    end     

    set(h2.fmin,'string',tmp.params.fmin); 
    set(h2.fmax,'string',tmp.params.fmax);
    set(h2.step,'string',tmp.params.step); 
    set(h2.smin,'string',tmp.params.smin); 
    set(h2.smax,'string',tmp.params.smax); 
    set(h2.ns,'string',tmp.params.ns);
    
%--- Instrument correction setting
    set(h2.poles,'string',tmp.params.poles);
    set(h2.data2,'string',tmp.params.data2);
    set(h2.save2,'string',tmp.params.save2);
          
case 'default3'
%-------------------------------- Preliminary analysis settings
%--- Spectral setting
    set(h2.w_spec,'string',tmp.params.w_spec); 
    set(h2.np,'string',tmp.params.np);
    set(h2.f_high_pass,'string',tmp.params.f_high_pass);
    set(h2.n_sta,'string',tmp.params.n_sta); 
    set(h2.spec_hour,'string',tmp.params.spec_hour);
%--- RMS setting
    set(h2.w_rms,'string',tmp.params.w_rms);
    set(h2.f1,'string',tmp.params.f1);
    set(h2.f2,'string',tmp.params.f2);
    set(h2.step_rms,'string',tmp.params.step_rms);
    set(h2.factor,'string',tmp.params.factor);  
%--- Polarization setting
    set(h2.w_pol,'string',tmp.params.w_pol);
    set(h2.f1_pol,'string',tmp.params.f1_pol);
    set(h2.f2_pol,'string',tmp.params.f2_pol);
    set(h2.step_pol,'string',tmp.params.step_pol);
    set(h2.factor_pol,'string',tmp.params.factor_pol);
%--- Detection setting
    set(h2.wshort,'string',tmp.params.wshort);
    set(h2.wlong,'string',tmp.params.wlong);
    set(h2.f1_det,'string',tmp.params.f1_det);
    set(h2.f2_det,'string',tmp.params.f2_det);
    set(h2.treshold,'string',tmp.params.treshold);
    set(h2.f1_down,'string',tmp.params.f1_down);
    set(h2.f2_down,'string',tmp.params.f2_down);
    set(h2.f1_up,'string',tmp.params.f1_up);
    set(h2.f2_up,'string',tmp.params.f2_up);
    
case 'default4'
%-------------------------------- Array analysis settings
%--- ZLC setting
    set(h2.w_zlc,'string',tmp.params.w_zlc);
    set(h2.f1_zlc,'string',tmp.params.f1_zlc);
    set(h2.f2_zlc,'string',tmp.params.f2_zlc);
    set(h2.step_zlc,'string',tmp.params.step_zlc);
    set(h2.vel_zlc,'string',tmp.params.vel_zlc);
    set(h2.splint,'string',tmp.params.splint);
    set(h2.maxlags,'string',tmp.params.maxlags);
%--- Radial Semblance/Semblance setting
    set(h2.w_sem,'string',tmp.params.w_sem);
    set(h2.f0,'string',tmp.params.f0);
    set(h2.Q,'string',tmp.params.Q);
    set(h2.vel_sem,'string',tmp.params.vel_sem);
    set(h2.n_sem,'string',tmp.params.n_sem);
    set(h2.step_sem,'string',tmp.params.step_sem);
    set(h2.topo_corr,'string',tmp.params.topo_corr);
    set(h2.freq_sem,'string',tmp.params.freq_sem);
    set(h2.xmin,'string',tmp.params.xmin);
    set(h2.xmax,'string',tmp.params.xmax);
    set(h2.ymin,'string',tmp.params.ymin);
    set(h2.ymax,'string',tmp.params.ymax);
    set(h2.zmin,'string',tmp.params.zmin);
    set(h2.zmax,'string',tmp.params.zmax);
%--- MUSIC setting
    set(h2.freq0,'string',tmp.params.freq0);
    set(h2.nfreq,'string',tmp.params.nfreq);
    set(h2.m,'string',tmp.params.m);
    set(h2.R,'string',tmp.params.R);
    set(h2.pmax,'string',tmp.params.pmax);
    set(h2.pinc,'string',tmp.params.pinc);
    set(h2.ninicio,'string',tmp.params.ninicio); 
    set(h2.lwin,'string',tmp.params.lwin);
    set(h2.advance,'string',tmp.params.advance); 
    set(h2.nsig,'string',tmp.params.nsig);
end
end