%%%%%%%%%%%%%Separate eigenvectors into signals and noise subspace%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N: number of stations
% T : length of the data window
% D : eigenvalues
% V : eigenvectors
% L : mofified eigenvalue, converted according to the method required
% E : corresponding eigenvectors
% flagsig: method to determine ranges of eigenvectors to use
% nfreq: the number of frequencies of the band of interest
% tprime: the analysis window for the MUSIC analysis (samples) 
% Nsig1, Nsig2: number of signal sources/ signal noise
% eigperc: the threshold value of the sources selection
% Nsig: the number of source computed  
% TMP_perc: the cumulative matrix/vector 

function [E,L,Nsig,TMP_perc,Nsig1,Nsig2]=signoise(row,row1,D,V,fo,flagsig,N,nfreq,tprime,Nsig1,Nsig2,eigperc)
%Order the eigenvalues and eigenvectors
[D,V]=Sort(N,D,V,row);
%Copy eigenarrays to complex array
for j=1:N
    for i=1:N
        E(i,j) = complex(V(i,j),V(i+N,j));
    end
    L(j) = D(j);
end
%Find the number of signals
Nsig=numsig(L,N,nfreq,tprime);
%Compute the cumulative matrix
TMP_cum=cumsum(L);
TMP_perc=TMP_cum./sum(L);
%Determine ranges of eigenvectros to use
if flagsig=='C'
    % CONVENTIONAL method;
    if Nsig1==0; Nsig1=1;end
    if Nsig2==0;Nsig2=N;end
elseif flagsig=='H'
    % HIGH RESOLUTION method
    if Nsig1==0;Nsig1=1;end
    if Nsig2==0;Nsig2=N;end
    for j=Nsig1:Nsig2
        L(j) = 1/D(j);
    end
elseif flagsig=='S'
    % MUSIC SIGNAL method
    if Nsig1==0;Nsig1=1;end
    if (Nsig2==0) | (Nsig2>=Nsig);Nsig2=Nsig;end;
    Nsig=Nsig2;
    for j=1:N
        L(j) = 1;
    end
elseif flagsig=='N'
    % AUTOMATIC method
    for j=1:N
        L(j) = 1;
    end
     delta=abs(TMP_perc-eigperc);
     [~,Nsig]=min(delta);
     Nsig1=Nsig+1;
     Nsig2=N+1;
end
end