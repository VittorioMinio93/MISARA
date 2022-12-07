%%%%%%%%%%%%%%Order the eigenvalues and eigenvectors%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N: number of stations
% NN: dimensions of eigenmatrixs
% D: eigenvalues matrix
% V: eigenvectors matrix
% s1: eigenvalues copy matrix
% V1: eigenvectors copy matrix

function [D,V]=Sort(N,D,V,row)
NN=2*N;
s1=D;
% Index the eigenvalues in descending order
% ix(i) stores the index of the i-th largest eigenvalue in D
for i=1:NN
    emax=-1000;
    for j=1:NN
        if s1(j)>emax
            emax=s1(j);
            jref=j;
        end 
    end
    ix(i)=jref;
    s1(jref)=-2000;    
end
% Swap if necessary
% eigenvalues are in pair
for i=1:2:NN
    ii=ix(i);
    jj=ix(i+1);
    if (V(1,ii)*V(N+1,jj))<0
        ix(i) = jj;
        ix(i+1) = ii;
    end
end
%Order the eigenvalues and eigenvectors
for j=1:N
    jref = ix(2*j-1);
    jref1= ix(2*j);
    s1(j) = D(jref);
    s1(j+N)= D(jref1);
    for i=1:NN
        V1(i,j) = V(i,jref);
        V1(i,j+N)=V(i,jref1);
    end
end
%Copy back to original array
for j=1:NN
    D(j) = s1(j);
    for i=1:NN
        V(i,j) = V1(i,j);      
    end
end
end