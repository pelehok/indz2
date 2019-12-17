function U_tU = du(tU,params)
n=length(tU);
U_tU=zeros(1,params.n);
step=params.step;
for i=1:n
    t=tU(i);
    nu=floor(t/step)+1;
    nu=min(nu,params.n);
    if (params.n==0)
        nu=1; 
    end
    t1=params.t0+(nu-1)*step;
    t2=t1+step;
    Ut=(t2-t)/step*params.U(nu)+(t-t1)/step*params.U(nu+1);
    U_tU(i)=Ut;
end 
end


