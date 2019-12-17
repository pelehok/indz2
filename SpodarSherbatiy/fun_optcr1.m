function [optcr] = fun_optcr1(u, params)
% Calculating the value of the optimization criteria 
params.U=u;

options = odeset('RelTol',1e-7,'AbsTol',1e-7);
soly=ode15s(@ode1_constr,[params.t0,params.Tend],params.y0,options,params);

optcr=trapz(soly.x,(soly.y(2,:)-params.yd(2)).^2);