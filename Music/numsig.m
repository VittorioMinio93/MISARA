%%%%%%%%%%%%%Finds the number of signals using the MDL method%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p,fp : number of stations (N)
% T,npts : length of window
% lambda : sorted eigenvalues (double precission)
% Nsig : number of signals, determined by MDL
% nfreq, fnfreq: the number of frequencies of the band of interest

function Nsig=numsig(lambda,p,nfreq,T)
npts=T;
mdlmin = 1.0e+20;
fp=p;
fnpts=npts;
fnfreq=nfreq;
for k=0:p-1
    sum = 0.0;
    sumlog = 0.0;
    fk=k;
    for i=k+1:p
        tmp=real(lambda(i));
        sum = sum + tmp;
        if tmp<0
            break
        end
        sumlog = sumlog + log(tmp);
    end
    amean = sum/(fp-fk);
    gmean = sumlog/(fp-fk);
    term=(fp-fk)*fnfreq*(log(amean)-gmean);
    mdl = term + 0.5*fk*(2.*fp-fk)*log(fnfreq);
    aic = term + fk*(2.*fp-fk);
    if mdl<mdlmin
        kmin=k;
        mdlmin = mdl;
    end   
end
Nsig = kmin;
sum  =0.0;
sum2 =0.0;
sum3 =0.0;
for i=1:p
    sum = lambda(i) + sum;
    sum2 = lambda(i)* lambda(i) + sum2;
    sum3 = lambda(i)*lambda(i)*lambda(i) + sum3;
end
beta1 = fp * sum2 - sum * sum;
beta1 = beta1 / ((fp-1.0) * sum * sum);
beta2 = sum3 - 1.5*sum*sum2+0.5*sum*sum*sum;
deno = (1.0/(fp*fp) - 1.5/fp + 0.5) * sum * sum * sum;
beta2 = 1.0 - beta2 / deno;
end