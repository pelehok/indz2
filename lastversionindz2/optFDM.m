function [cr, grad] = optFDM(u, params)
params.U=u;
options = odeset('RelTol',1e-7,'AbsTol',1e-7);
soly=ode15s(@ode,[params.t0,params.T],params.y0,options,params);
cr=trapz(soly.x,(soly.y(2,:)-params.yd(2)).^2);
grad = FDM(u,params);

