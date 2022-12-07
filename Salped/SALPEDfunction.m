%%%%%%%Salped algorithm (Exclusively Long Period signals detection)%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s: seismic trace
% fs: sampling interval
% f1: the minimum frequency for the central band of interest(Hz)
% f2: the maximum frequency for the central band of interest(Hz)
% f1_down: the minimum frequency for the lower band of no interest(Hz)
% f2_down: the maximum frequency for the lower band of no interest(Hz)
% f1_up: the minimum frequency for the upper band of no interest(Hz)
% f2_up: the maximum frequency for the upper band of no interest(Hz)
% CF:the Characteristic function
% tCF: time vector of the Characteristic function

function [CF,tCF]=SALPEDfunction(s,fs,f1,f2,f1_down,f2_down,f1_up,f2_up)
%Rescale
s=s-mean(s);
T_fir=4;%Hamming window (s)
longitud_fir=ceil(T_fir*fs);%Hamming window(samples)

%Lower subband 
BL_i = fir1(longitud_fir,[f1_down f2_down]/(fs/2));
BL_q=imag(hilbert(BL_i));
%Central subband 
BC_i = fir1(longitud_fir,[f1 f2]/(fs/2)); 
BC_q=imag(hilbert(BC_i));
%Upper subband 
BU_i = fir1(longitud_fir,[f1_up f2_up]/(fs/2));  
BU_q=imag(hilbert(BU_i));

%Band-pass filtering and extracting of the spindle-shaped envelopes for the Lower,Central and Upper subbands 
xL=filter(BL_i,1,s); xL_q=filter(BL_q,1,s); eL=sqrt(xL.^2 + xL_q.^2);
xC=filter(BC_i,1,s); xC_q=filter(BC_q,1,s); eC=sqrt(xC.^2 + xC_q.^2);
xU=filter(BU_i,1,s); xU_q=filter(BU_q,1,s); eU=sqrt(xU.^2 + xU_q.^2);

B_decima=fir1(longitud_fir,0.33/fs/2); % Frequency filtering: 0.33Hz
Factor_dec = round(fs/2);              % To decimate the envelopes: 2 samples per sec

%Low-pass fltering and decimation of the envelopes
eL=filter(B_decima,1,eL); eL=decimate(eL,Factor_dec);
eC=filter(B_decima,1,eC); eC=decimate(eC,Factor_dec);
eU=filter(B_decima,1,eU); eU=decimate(eU,Factor_dec);

% Provide robustness against noise.
B_lp = fir1(400,0.005);
e_lp = filtfilt(B_lp,1,eL+eC+eU);%Low-pass filtering (0.005 Hz)
%Normalization
eLp=eL./e_lp;
eCp=eC./e_lp;
eUp=eU./e_lp;

% Application of threshold to obtain the noise-robust envelopes
THR = 0.4;
zL=eLp.*(eLp>THR)+THR.*(eLp<=THR);
zC=eCp.*(eCp>THR)+THR.*(eCp<=THR);
zU=eUp.*(eUp>THR)+THR.*(eUp<=THR);


% Discriminant detector 
fs_reduc=2;                                     % 2 samples per sec
B_hat=mexihat(-4,4,15.5*fs_reduc);              % filter to detect LPs (LPs duration expected : 4 sec)
B_hat=B_hat.*(B_hat>0) + 1.6*B_hat.*(B_hat<=0); % Events with longer durations are punished with a negative impulse response (penalty function: 1.6)

% Detection with  "Mexican Hat" filter
yL = filter(B_hat,1,zL);
yC = filter(B_hat,1,zC);
yU = filter(B_hat,1,zU);

% Application of threshold 
THR=0;
yL=yL.*(yL>THR)+THR.*(yL<=THR);
yC=yC.*(yC>THR)+THR.*(yC<=THR);
yU=yU.*(yU>THR)+THR.*(yU<=THR);

% Computing of the temporal vector
tCF =(0:(length(yC)-1))/fs_reduc - 11.75;      % compensate for time delay  (2+2+7.75 sec)  

% Penlize false positive detections with spectral penalty factors
M=2;    
N=4;  
%Computing of the characterisic function 
CF=yC-M*yL-N*yU;

CF(CF<0)=0;% Replace the negative impulses of the characterisic function with zero values  
CF(tCF<5)=0; % Replace the unreliable part of the characterisic function with zero values  
end