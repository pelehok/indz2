function [ c,constr, grad] = constrPsi1AM( u,params )
params.U=u;
options = odeset('RelTol',1e-7,'AbsTol',1e-7);
soly=ode15s(@ode1_constr,[params.t0,params.Tend],params.y0,options,params);
c=[];
constr = trapz(soly.x',(soly.y(1,:)-params.yy(1)+abs(soly.y(1,:)-params.yy(1))).^2);

% syms y1 yy
% grad = gradient((y1-yy+abs(y1-yy))^2, [y1]);

count = 100;
dfdb=0*ones(2,count);
grad = 0*ones(1,params.n+1);
for k=params.n+1:-1:1
    params.k=k;
    linsp = linspace (params.Tend, params.t0, count);
    solmiu=ode15s(@odeMiuConstr,linsp,params.miu0,options,params);  
    miu=deval(solmiu, linsp);
    yy=deval(params.sol, linsp);
    index=count;
    for q=params.Tend:-(params.Tend-params.t0)/(count-1):params.t0
        dfdb([1, 2], index)=find_dfdb(q,params,yy(:,index));
        index=index-1;
    end    
    grad(params.n-k+2) = -trapz(linsp', dfdb(1,:).*miu(1,:)+dfdb(2,:).*miu(2,:));
end
end

