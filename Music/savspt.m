%%%%%%%%%%%%%Return the real and imaginary parts of the spectrums %%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CF: index of frequency analysis
% M: Hamming value
% FLAGLA: Active/Deactive the normalization of the spectrums 
% cu: spectrums matrix
% AMP, PHI: real and imaginary parts of the spectrums
% N: number of stations
function [PHI,AMP]=savspt(CF,M,FLAGLA,N,cu)
PHI=zeros(1,N);%Initialize imaginary part of the spectrum
AMP=zeros(1,N);%Initialize real part of the spectrum
pi=4*tan(1);
%Loop through the number of stations
for ista=1:N
    k=1;
    %Loop for each frequency in frequency smoothing band
    for j=CF-M:CF+M
        %Active/deactive normalization
        if (FLAGLA=='Y') | (FLAGLA=='y')
            denom=abs(cu(j,ista));
            %Compute real and imaginary parts of the spectrums
            AMP(k,ista)=real(cu(j,ista))/denom;
            PHI(k,ista)=imag(cu(j,ista))/denom; 
        else
            %Compute real and imaginary parts of the spectrums
            AMP(k,ista)=real(cu(j,ista));
            PHI(k,ista)=imag(cu(j,ista));
        end
        k=k+1;
    end
end
end