%%%%%%%%%%%%% Control and reset automatically user parameters in the Home window%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio

%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tmp: user parameters
%h2: UI control parameters
%%%%%%%%%%%h2.variables/params.variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------- General settings
%--- General info
% sta_ref: the name of the reference station 
%-------------------------------- Instrumental settings
%--- Coordinate system/Beam Pattern setting
% fmin: the minimum frequency for the Beam pattern analysis (Hz)
% fmax: the maximum frequency for the Beam pattern analysis(Hz)
% step: the frequency step for the Beam pattern analysis (Hz)
% smin: the minimum slowness for the Beam pattern analysis (s/km)
% smax: the miximum slowness for the Beam pattern analysis (s/km)
% ns: the dimensions of the square grid slowness for the Beam pattern analysis
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

function setdef(hobj,~)
%% Declaration global parameters
global h2
%% load user parameters
tmp= load('user_values.mat','params');

%% Control  and set parameters
%%Execute one of four groups of settings
switch get(hobj,'Tag')
case 'setdef1'
%-------------------------------- General settings
%--- General info
sta_ref=deblank(get(h2.sta_ref,'string'));
if (size(sta_ref,2)==0);set(h2.sta_ref,'string',tmp.params.sta_ref);end

case 'setdef2' 
%-------------------------------- Instrumental settings
%--- Coordinate system/Beam Pattern setting
fmin=str2num(get(h2.fmin,'string'));
fmax=str2num(get(h2.fmax,'string')); 
step=str2num(get(h2.step,'string'));
smin=str2num(get(h2.smin,'string'));
smax=str2num(get(h2.smax,'string'));
ns=str2num(get(h2.ns,'string'));
if (size(fmin,1)==0);set(h2.fmin,'string',tmp.params.fmin);end
if (size(fmax,1)==0);set(h2.fmax,'string',tmp.params.fmax);end
if (size(step,1)==0);set(h2.step,'string',tmp.params.step);end 
if (size(smin,1)==0);set(h2.smin,'string',tmp.params.smin);end
if (size(smax,1)==0);set(h2.smax,'string',tmp.params.smax);end 
if (size(ns,1)==0);set(h2.ns,'string',tmp.params.ns);end

case 'setdef3' 
%-------------------------------- Preliminary analysis settings
%--- Spectral setting
w_spec=str2num(get(h2.w_spec,'string'));
np=str2num(get(h2.np,'string'));
f_high_pass=str2num(get(h2.f_high_pass,'string'));
n_sta=str2num(get(h2.n_sta,'string'));
spec_hour=str2num(get(h2.spec_hour,'string'));
if (size(w_spec,1)==0);set(h2.w_spec,'string',tmp.params.w_spec);end
if (size(np,1)==0);set(h2.np,'string',tmp.params.np);end
if (size(f_high_pass,1)==0);set(h2.f_high_pass,'string',tmp.params.f_high_pass);end
if (size(n_sta,1)==0);set(h2.n_sta,'string',tmp.params.n_sta);end
if (size(spec_hour,1)==0);set(h2.spec_hour,'string',tmp.params.spec_hour);end

%--- RMS setting
w_rms=str2num(get(h2.w_rms,'string'));
f1=str2num(get(h2.f1,'string'));
f2=str2num(get(h2.f2,'string'));
step_rms=str2num(get(h2.step_rms,'string'));
factor=str2num(get(h2.factor,'string'));
if (size(w_rms,1)==0);set(h2.w_rms,'string',tmp.params.w_rms);end
if (size(f1,1)==0);set(h2.f1,'string',tmp.params.f1);end
if (size(f2,1)==0);set(h2.f2,'string',tmp.params.f2);end
if (size(step_rms,1)==0);set(h2.step_rms,'string',tmp.params.step_rms);end
if (size(factor,1)==0);set(h2.factor,'string',tmp.params.factor);end

%--- Polarization setting
w_pol=str2num(get(h2.w_pol,'string'));
f1_pol=str2num(get(h2.f1_pol,'string'));
f2_pol=str2num(get(h2.f2_pol,'string'));
step_pol=str2num(get(h2.step_pol,'string'));
factor_pol=str2num(get(h2.factor_pol,'string'));
if (size(w_pol,1)==0);set(h2.w_pol,'string',tmp.params.w_pol);end
if (size(f1_pol,1)==0);set(h2.f1_pol,'string',tmp.params.f1_pol);end
if (size(f2_pol,1)==0);set(h2.f2_pol,'string',tmp.params.f2_pol);end
if (size(step_pol,1)==0);set(h2.step_pol,'string',tmp.params.step_pol);end
if (size(factor_pol,1)==0);set(h2.factor_pol,'string',tmp.params.factor_pol);end

%--- Detection setting
wshort=str2num(get(h2.wshort,'string'));
wlong=str2num(get(h2.wlong,'string'));
f1_det=str2num(get(h2.f1_det,'string'));
f2_det=str2num(get(h2.f2_det,'string'));
treshold=str2num(get(h2.treshold,'string'));
f1_down=str2num(get(h2.f1_down,'string'));
f2_down=str2num(get(h2.f2_down,'string'));
f1_up=str2num(get(h2.f1_up,'string'));
f2_up=str2num(get(h2.f2_up,'string'));
if (size(wshort,1)==0);set(h2.wshort,'string',tmp.params.wshort);end
if (size(wlong,1)==0);set(h2.wlong,'string',tmp.params.wlong);end
if (size(f1_det,1)==0);set(h2.f1_det,'string',tmp.params.f1_det);end
if (size(f2_det,1)==0);set(h2.f2_det,'string',tmp.params.f2_det);end
if (size(treshold,1)==0);set(h2.treshold,'string',tmp.params.treshold);end
if (size(f1_down,1)==0);set(h2.f1_down,'string',tmp.params.f1_down);end
if (size(f2_down,1)==0);set(h2.f2_down,'string',tmp.params.f2_down);end
if (size(f1_up,1)==0);set(h2.f1_up,'string',tmp.params.f1_up);end
if (size(f2_up,1)==0);set(h2.f2_up,'string',tmp.params.f2_up);end


case 'setdef4'
%-------------------------------- Array analysis settings
%--- ZLC setting
w_zlc=str2num(get(h2.w_zlc,'string'));
f1_zlc=str2num(get(h2.f1_zlc,'string'));
f2_zlc=str2num(get(h2.f2_zlc,'string'));
step_zlc=str2num(get(h2.step_zlc,'string'));
vel_zlc=str2num(get(h2.vel_zlc,'string'));
splint=str2num(get(h2.splint,'string'));
maxlags=str2num(get(h2.maxlags,'string'));
if (size(w_zlc,1)==0);set(h2.w_zlc,'string',tmp.params.w_zlc);end
if (size(f1_zlc,1)==0);set(h2.f1_zlc,'string',tmp.params.f1_zlc);end
if (size(f2_zlc,1)==0);set(h2.f2_zlc,'string',tmp.params.f2_zlc);end
if (size(step_zlc,1)==0);set(h2.step_zlc,'string',tmp.params.step_zlc);end
if (size(vel_zlc,1)==0);set(h2.vel_zlc,'string',tmp.params.vel_zlc);end 
if (size(splint,1)==0);set(h2.splint,'string',tmp.params.splint);end
if (size(maxlags,1)==0);set(h2.maxlags,'string',tmp.params.maxlags);end 

%--- Radial Semblance/Semblance setting
w_sem=str2num(get(h2.w_sem,'string'));
f0=str2num(get(h2.f0,'string'));
Q=str2num(get(h2.Q,'string'));
vel_sem=str2num(get(h2.vel_sem,'string'));
n_sem=str2num(get(h2.n_sem,'string'));
step_sem=str2num(get(h2.step_sem,'string'));
topo_corr=deblank(get(h2.topo_corr,'string'));
freq_sem=deblank(get(h2.freq_sem,'string'));
xmin=str2num(get(h2.xmin,'string'));
xmax=str2num(get(h2.xmax,'string'));
ymin=str2num(get(h2.ymin,'string'));
ymax=str2num(get(h2.ymax,'string'));
zmin=str2num(get(h2.zmin,'string'));
zmax=str2num(get(h2.zmax,'string'));
if (size(w_sem,1)==0);set(h2.w_sem,'string',tmp.params.w_sem);end
if (size(f0,1)==0);set(h2.f0,'string',tmp.params.f0);end 
if (size(Q,1)==0);set(h2.Q,'string',tmp.params.Q);end
if (size(vel_sem,1)==0);set(h2.vel_sem,'string',tmp.params.vel_sem);end
if (size(n_sem,1)==0);set(h2.n_sem,'string',tmp.params.n_sem);end 
if (size(step_sem,1)==0);set(h2.step_sem,'string',tmp.params.step_sem);end
if (size(topo_corr,1)==0);set(h2.topo_corr,'string',tmp.params.topo_corr);end
if (size(freq_sem,2)==0);set(h2.freq_sem,'string',tmp.params.freq_sem);end
if (size(xmin,1)==0);set(h2.xmin,'string',tmp.params.xmin);end
if (size(xmax,1)==0);set(h2.xmax,'string',tmp.params.xmax);end
if (size(ymin,1)==0);set(h2.ymin,'string',tmp.params.ymin);end 
if (size(ymax,1)==0);set(h2.ymax,'string',tmp.params.ymax);end
if (size(zmin,1)==0);set(h2.zmin,'string',tmp.params.zmin);end
if (size(zmax,1)==0);set(h2.zmax,'string',tmp.params.zmax);end

%--- MUSIC setting
freq0=str2num(get(h2.freq0,'string'));
nfreq=str2num(get(h2.nfreq,'string'));
m=str2num(get(h2.m,'string'));
R=str2num(get(h2.R,'string'));
pmax=str2num(get(h2.pmax,'string'));
pinc=str2num(get(h2.pinc,'string'));
ninicio=str2num(get(h2.ninicio,'string'));
lwin=str2num(get(h2.lwin,'string'));
advance=str2num(get(h2.advance,'string'));
nsig=str2num(get(h2.nsig,'string'));
if (size(freq0,1)==0);set(h2.freq0,'string',tmp.params.freq0);end  
if (size(nfreq,1)==0);set(h2.nfreq,'string',tmp.params.nfreq);end 
if (size(m,1)==0);set(h2.m,'string',tmp.params.m);end
if (size(R,1)==0);set(h2.R,'string',tmp.params.R);end
if (size(pmax,1)==0);set(h2.pmax,'string',tmp.params.pmax);end  
if (size(pinc,1)==0);set(h2.pinc,'string',tmp.params.pinc);end 
if (size(ninicio,1)==0);set(h2.ninicio,'string',tmp.params.ninicio);end
if (size(lwin,1)==0);set(h2.lwin,'string',tmp.params.lwin);end
if (size(advance,1)==0);set(h2.advance,'string',tmp.params.advance);end  
if (size(nsig,1)==0);set(h2.nsig,'string',tmp.params.nsig);end  
end
end