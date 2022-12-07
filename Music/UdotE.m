%%%%%%%%%%%%%%Compute the directional function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N: number of stations
% kx,ky: wave numbers
% E: eigenvectors selected matrix
% L: projection length
% x: station coordinates matrix
% sigflag: method to determine ranges of eigenvectors to use
% Nsig: the number of source computed 
% Nsig1, Nsig2: number of signal sources/ signal noise

function UE=UdotE(N,kx,ky,E,L,x,row1,sigflag,Nsig,Nsig1,Nsig2)
% This subroutine computes the dot product in the directional function.
% The directional vector used in this program is: exp{-i(kx*x+ky*y)}/sqrt(N), N = no. of stations
% Compute  {Hermitain(UE) * (E)}
% UE is the directional vector: exp{-(kx*x+ky*y)}
% E is the eigenv matrix
sum1=0;
if (sigflag=='N') | (sigflag=='S')
    for j=1:Nsig1-1
        sum=complex(0,0);
        for i=1:N
            arg = kx*x(1,i) + ky*x(2,i);
	        carg =complex(0,arg);
            sum = sum + exp(carg)*E(i,j);
        end
        sum1 = sum1 + L(j)*(abs(sum)^2); 
    end
    UE = 1.0 - sum1/N;
else
    for j=Nsig1:Nsig2
        sum=complex(0,0);
        for i=1:N
            arg = kx*x(1,i) + ky*x(2,i);
	        carg =complex(0,arg);
            sum = sum + exp(carg)*E(i,j);
        end
        sum1 = sum1 + L(j)*(abs(sum)^2);  
    end
    UE = sum1/N;   
end
end