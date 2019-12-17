function dy = ode1_constr(t,y,params)

if params.n==0 nu=1; else
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

dy=[params.p(1)*y(1)+params.p(4)*y(1)*(1-exp(-params.p(5)*y(1)))*y(2)/...
    (params.p(6)+params.p(7)*y(2));...
    params.p(2)*y(2)+params.p(3)*y(2)*y(2)+params.p(8)*params.p(4)*y(1)*...
    (1-exp(-params.p(5)*y(1)))*y(2)/(params.p(6)+params.p(7)*y(2))]; 

end