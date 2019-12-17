function dy = ode(t,y,params)
if params.n==0
    countofu=1; 
else
step=params.step;
countofu=floor(t/step)+1;
countofu=min(countofu,params.n);
end
        t1=params.t0+(countofu-1)*step;
        t2=t1+step;
        Ut=(t2-t)/step*params.U(countofu)+(t-t1)/step*params.U(countofu+1); 

params.p(params.j)=Ut;
   dy=zeros(2,1);
   dy(1)=params.p(1)*y(1)+params.p(4)*y(1)*(1-exp(-params.p(5)*y(1)))*y(2)/...
    (params.p(6)+params.p(7)*y(2));
dy(2)=params.p(2)*y(2)+params.p(3)*y(2)*y(2)+params.p(8)*params.p(4)*y(1)*...
    (1-exp(-params.p(5)*y(1)))*y(2)/(params.p(6)+params.p(7)*y(2));
end

