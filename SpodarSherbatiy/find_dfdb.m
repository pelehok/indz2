function [dfdb] = find_dfdb(t,params, yi)

if params.n==0 
    nu=1; else
lts=params.lts;
nu=floor(t/lts)+1;
nu=min(nu,params.n);
end
switch params.typeU
    case 'const' 
       Ut=params.U(nu);
    case 'linear'
        t1=params.t0+(nu-1)*lts;
        t2=t1+lts;
        Ut=(t2-t)/lts*params.U(nu)+(t-t1)/lts*params.U(nu+1);
 
    otherwise
        error('Unexpected approximation type')
end
params.p(params.j)=Ut;

if (params.n==0) 
    coef=1; 
else
    if (nu==params.k) 
        coef = (t2-t)/lts;       
    %elseif(params.k==params.n+1) 
    elseif (nu+1==params.k && nu ~= t/lts+1)
        coef = (t-t1)/lts;         
    else
    coef=0;
    end
end

dfdb=[-(yi(2).*params.p(4).*yi(1).*(1-exp(-params.p(5)*yi(1))))./(params.p(6)+...
    params.p(7)*yi(2))./(params.p(6)+params.p(7).*yi(2)).*coef;
    -(yi(2).*params.p(8)*params.p(4).*yi(1).*(1-exp(-params.p(5).*yi(1))))./(params.p(6)+...
    params.p(7)*yi(2))./(params.p(6)+params.p(7)*yi(2)).*coef];

dfdb=[(yi(1).*(1-exp(-params.p(5).*yi(1)))*...
    yi(2)./(params.p(6)+params.p(7)*yi(2))).*coef;
    (params.p(8)*yi(1).*(1-exp(-params.p(5).*yi(1)))*...
    yi(2)./(params.p(6)+params.p(7)*yi(2))).*coef];

end
