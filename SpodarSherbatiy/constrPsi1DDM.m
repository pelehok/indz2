function [ c,constr, grad] = constrPsi1DDM( u,params )
params.U=u;
options = odeset('RelTol',1e-7,'AbsTol',1e-7);
soly=ode15s(@ode1_constr,[params.t0,params.Tend],params.y0,options,params);
c=[];
constr = trapz(soly.x',(soly.y(1,:)-params.yy(1)+abs(soly.y(1,:)-params.yy(1))).^2);

% syms y1 yy
% grad = gradient((y1-yy+abs(y1-yy))^2, [y1]);

grad = 0*ones(1,params.n+1);
for k=1:params.n+1
    params.k=k;
    solz=ode15s(@odeZ,[params.t0, params.Tend],params.z0,options,params);
    linsp = linspace (params.t0, params.Tend, 100);
    zi=deval(solz, linsp);
    yy=deval(soly, linsp);
    dg1dy=[2*(abs(yy(1,:)-params.yy(1))+yy(1,:)-params.yy(1))'.*...
    (sign(yy(1,:)-params.yy(1))+1)', (0*ones(1,size(yy,2)))'];
    grad(k) = trapz(linsp', dg1dy(:,1).*zi(1,:)');
end
end

