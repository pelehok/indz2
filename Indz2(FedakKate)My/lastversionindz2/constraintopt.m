function [c,constraint]= constraintopt( u,params )
params.U=u;
options = odeset('RelTol',1e-7,'AbsTol',1e-7);
[t,y]=ode15s(@ode,[params.t0,params.T],params.y0,options,params);
c=[];
constraint = trapz(t,(y(:,2)-params.yMax(2)+abs(y(:,2)-params.yMax(2))).^2);
end

