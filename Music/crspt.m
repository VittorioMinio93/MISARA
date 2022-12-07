%%%%%%%%%%%%%Estimation of the spatial cross-spectral matrix %%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N: number of stations
% M:Hamming value
% phi, amp: real and imaginary parts of the spectrums
% a: Hamming points
% tt: the focusing matrix
% temp: the covariance matrix
% rr: the estimate of the cross-spectral matrix in partitioned tilted format
%       | real      |-imaginary |
%       |-----------------------| = rr
%       | imaginary |   real    |
% re: real part
% im:  imaginary part

function rr=crspt (N,M,phi,amp,a,tt)
temp=zeros(2*N,N);%Initialize the covariance matrix  
rr=complex(zeros(N,N));% Initialize the cross-spectral matrix
% Loop through all stations pairs
for i=1:N
    re=0;
    %Diagonal elements (autospectrum; real number)
    %smoothing over 2*M+1 frequencies
    for l=1:2*M+1
        re = (amp(l,i)*amp(l,i)+phi(l,i)*phi(l,i))*a(l) + re;
    end
    temp(i,i) = re;
    temp(i+N,i) = 0.0;
    %Off-diagonal elements (crossspectrum; complex number)
    for j=i+1:N
         re = 0.0;
         im = 0.0;
         for l=1:2*M+1
             %Form Cross-spectral matrix estimate : dft(i) * Hermitian (dft(j))
             re =(  amp(l,i)*amp(l,j)+phi(l,i)*phi(l,j)  ) * a(l)+ re;
             im =(  amp(l,j)*phi(l,i)-phi(l,j)*amp(l,i)  ) * a(l)+ im;
         end
         temp(i,j) = complex(re,im);
         temp(j,i) = conj(temp(i,j));  
    end
end
% Apply focusing transformation: tt*rr*Herm(tt)
% The transformation matrix ttis NxN.
% tt gives approximate transformation from frequency f to frequency fo
for i=1:N
    for j=1:N
        temp1=complex(0);
        for k=1:N
            for l=1:N
                temp1 = temp1 + tt(i,k)*conj(tt(j,l))*temp(k,l);
            end
        end
        rr(i,j)=real(temp1);
        rr(i+N,j)=imag(temp1);
        rr(i+N,j+N)=rr(i,j);
        rr(i,j+N)=-rr(i+N,j);
    end
end
end