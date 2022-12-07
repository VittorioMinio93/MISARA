%%%%%%%%%%%%%%Find the main peaks in the fk-power spectrum%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ns: dimensions of slowness grid
% ds: the size step of the grid slowness (s/km) 
% sx0, sy0: the limits of the grid slowness (s/km)
% Nsig: the number of sources computed 
% d: power spectrum
% peak: main peaks of the fk-power spectrum

function peak=fpeak(ns,ds,sx0,sy0,Nsig,d)
peak=zeros(4,Nsig);%Initialize the main peak matrix
% Loop over the slwoness grid except over the edges
for i=2:ns-1
    for j=2:ns-1
        iflag=0;
        iis = i - 1;
        iil = i + 1;
        jjs = j - 1;
        jjl = j + 1;
        if i==1;iis = 1;end
        if i==ns; iil = ns;end
        if j==1; jjs = 1;end
        if j==ns; jjl = ns;end
        for ii=iis:iil
            for jj=jjs:jjl
                if d(ii,jj)>d(i,j)
                    iflag=1;
                end
            end
        end
        if iflag==0
            kref = Nsig+1;
            for kk=Nsig:-1:1
                if d(i,j)>peak(1,kk)
                    kref=kk;
                end
            end
            if kref<=Nsig
                if kref<Nsig
                    for kk=Nsig-1:-1:kref
                        for ll=1:4
                            peak(ll,kk+1) = peak(ll,kk);
                        end
                    end
                end
                sx = (i-1)*ds + sx0;
                sy = (j-1)*ds + sy0;
                peak(1,kref) = d(i,j);
                peak(2,kref) = sx;
                peak(3,kref) = sy;
                test=sx*sx+sy*sy;
                if test>0
                    peak(4,kref) = 1./sqrt(sx*sx+sy*sy);
                else
                    peak(4,kref) = 1.0e+20;
                end
            end
        end   
    end   
end
end