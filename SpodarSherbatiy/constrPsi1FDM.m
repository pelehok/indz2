function [ c,constr, grad] = constrPsi1FDM( u,params )
params.U=u;
options = odeset('RelTol',1e-7,'AbsTol',1e-7);
soly=ode15s(@ode1_constr,[params.t0,params.Tend],params.y0,options,params);
c=[];
constr = trapz(soly.x',(soly.y(1,:)-params.yy(1)+abs(soly.y(1,:)-params.yy(1))).^2);
[ccc,grad] = Differences (u,params);

