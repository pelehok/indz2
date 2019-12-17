function U_tU = Values_U_tU(tU,params)
%Exersize 2.1.1 - right hand side of equation 
n=length(tU);
U_tU=zeros(1,n);
lts=params.lts;
for i=1:n
    t=tU(i);
    nu=floor(t/lts)+1;
    nu=min(nu,params.n);
    if (params.n==0) nu=1; end
    switch params.typeU
        case 'const' 
           Ut=params.U(nu);
        case 'linear'
            t1=params.t0+(nu-1)*lts;
            t2=t1+lts;
            Ut=(t2-t)/lts*params.U(nu)+(t-t1)/lts*params.U(nu+1);
        otherwise
            error('Unexpected opproximation type')
    end %switch
    U_tU(i)=Ut;
end %for

end
