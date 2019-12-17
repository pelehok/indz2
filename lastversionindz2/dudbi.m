function dudb = dudbi(t,params,y)

if params.n==0 
    nu=1; else
    step=params.step;
    nu=floor(t/step)+1;
    nu=min(nu,params.n);
end

t1=params.t0+(nu-1)*step;
t2=t1+step;
Ut=(t2-t)/step*params.U(nu)+(t-t1)/step*params.U(nu+1);
params.p(params.j)=Ut;

if (params.n==0) 
    coef=1;
else
    if (nu==params.k) 
        coef = (t2-t)/step;
    elseif (nu+1==params.k && nu ~= t/step+1)
        coef = (t-t1)/step;
    else
    coef=0;
    end
end

dudb=[(params.p(4).*(y(1)^2).*exp(-params.p(5).*y(1))*...
    y(2)./(params.p(6)+params.p(7)*y(2))).*coef;
    (params.p(8)*params.p(4).*(y(1)^2).*exp(-params.p(5).*y(1))*...
    y(2)./(params.p(6)+params.p(7)*y(2))).*coef];
end

