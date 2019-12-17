function [criteria] = criteria(u, params)
% Calculating the value of the optimization criteria 
params.U=u;
options = odeset('RelTol',1e-7,'AbsTol',1e-7);
soly=ode15s(@ode,[params.t0,params.T],params.y0,options,params);
criteria=trapz(soly.x,(soly.y(1,:)-params.yd(1)).^2);
end

