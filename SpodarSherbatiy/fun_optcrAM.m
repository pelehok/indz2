function [optcr, grad] = fun_optcrAM(u, params)
% Calculating the value of the optimization criteria 
params.U=u;

options = odeset('RelTol',1e-7,'AbsTol',1e-7);
soly=ode15s(@ode1_constr,[params.t0,params.Tend],params.y0,options,params);

optcr=trapz(soly.x,(soly.y(2,:)-params.yd(2)).^2);

count = 100;
dfdb=0*ones(2,count);
grad = 0*ones(1,params.n+1);
linsp = linspace (params.Tend, params.t0, count);
solmiu=ode15s(@odeMiu,linsp,params.miu0,options,params);  
miu=deval(solmiu, linsp);
yy=deval(soly, linsp);
for k=params.n+1:-1:1
    params.k=k;  
    index=count;
    for q=params.Tend:-(params.Tend-params.t0)/(count-1):params.t0
        dfdb([1, 2], index)=find_dfdb(q,params,yy(:,index));
        index=index-1;
    end    
    grad(params.n-k+2) = -trapz(linsp', dfdb(1,:).*miu(1,:)+dfdb(2,:).*miu(2,:));
end
