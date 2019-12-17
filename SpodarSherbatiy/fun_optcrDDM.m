function [optcr, grad] = fun_optcrDDM(u, params)
% Calculating the value of the optimization criteria 
params.U=u;

options = odeset('RelTol',1e-7,'AbsTol',1e-7);
soly=ode15s(@ode1_constr,[params.t0,params.Tend],params.y0,options,params);

optcr=trapz(soly.x,(soly.y(2,:)-params.yd(2)).^2);

dg1dy=Find_dg1dy (soly.y);
grad = 0*ones(1,params.n+1);
linsp = linspace (params.t0, params.Tend, 100);
yy=deval(soly, linsp);
dg1dy=Find_dg1dy(yy);
for k=1:params.n+1
    params.k=k;
    solz=ode15s(@odeZ,[params.t0, params.Tend],params.z0,options,params); 
    zi=deval(solz, linsp);       
    grad(k) = trapz(linsp', dg1dy(:,2).*zi(2,:)');
end
