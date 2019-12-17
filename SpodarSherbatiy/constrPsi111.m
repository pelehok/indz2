function [ constr ] = constrPsi111( u,params )
params.U=u;
options = odeset('RelTol',1e-7,'AbsTol',1e-7);
[t,y]=ode15s(@ode1_constr,[params.t0,params.Tend],params.y0,options,params);
constr = trapz(t,(y(:,1)-params.yy(1)+abs(y(:,1)-params.yy(1))).^2);
end

