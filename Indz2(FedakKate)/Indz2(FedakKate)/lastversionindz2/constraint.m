function [constraint]  = constraint( u,params )
params.U=u;
options = odeset('RelTol',1e-7,'AbsTol',1e-7);
[t,y]=ode15s(@ode,[params.t0,params.T],params.y0,options,params);
constraint = trapz(t,(y(:,1)-params.yMax(1)+abs(y(:,1)-params.yMax(1))).^2);
end

