%%%%%%%%%%%%%  Set the Hamming window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%This routine constructs a Hamming window array "a" of 2M+1 points.
%If flag =1 then exclude central 2M+1 frequencies (BR method)
%%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,suma2]=haming(M)
         suma=0;
         if M==0
           a=1;
           suma2=1;
         else
             for i=1:M
                 a(M+1-i)=0.54 + 0.46*cos(pi*(i/M));
                 a(i+M+1) = a(M+1-i);
                 suma = suma + 2.*a(M+1-i);
             end
             a(M+1) = 1;
             suma = suma + a(M+1);
             suma2 = 0.0;
             for i=1:2*M+1
                 a(i) = a(i)/suma;
                 suma2 = suma2 + a(i)*a(i);  
             end
             
         end
end