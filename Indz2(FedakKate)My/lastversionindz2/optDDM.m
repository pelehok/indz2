function [optcr, constreint] = optDDM(u, params)
    params.U=u;
    options = odeset('RelTol',1e-7,'AbsTol',1e-7);
    soly=ode15s(@ode,[params.t0,params.T],params.y0,options,params);
    optcr=trapz(soly.x,(soly.y(1,:)-params.yd(1)).^2);
    constreint = zeros(1,params.n+1);
    linsp = linspace (params.t0, params.T, 100);
    y=deval(soly, linsp);
    g0=[(zeros(1,size(y,2)))', 2*y(2,:)'-1];
    for j=1:params.n+1
        params.k=j;
        solDDM=ode15s(@odeDDM,[params.t0,params.T],[0 0],options,params);
        zi=deval(solDDM, linsp);  
        constreint(j) = trapz(linsp',g0(:,2).*zi(2,:)');
    end
end