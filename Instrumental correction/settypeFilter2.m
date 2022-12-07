%%%%%%%%%%%%% Filters of the seismic traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% z: amplitude of the seismic trace
% delim: type of filter
% f1, f2: frequencies filters
% fs: sampling rate 

function z=settypeFilter2(z,delim,f1,f2,fs)
         z=z;
         if delim==1
             %Low-pass filter
             fn=fs/2;
             W=f2/fn;
             [x y]=butter(2,W,'low');
             z=filter(x,y,z);
         elseif delim==2
             %High-pass filter 
             fn=fs/2;
             W=f1/fn;
            [x y]=butter(2,W,'high');
            z=filter(x,y,z);
         elseif delim==3
             %Band-pass filter 
             fn=fs/2;
             W=[f1 f2]/fn;
             [x,y]=butter(2,W);
             z=filter(x,y,z);
         end
%          z=num2cell(z,[1,size(z,2)]);
end