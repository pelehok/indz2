function [cr, constreint] = optAM(u, params)
    params.U=u;
    options = odeset('RelTol',1e-7,'AbsTol',1e-7);
    soly=ode15s(@ode,[params.t0,params.T],params.y0,options,params);
    cr=trapz(soly.x,(soly.y(1,:)-params.yd(1)).^2);
    count = 100;
    dfdb=zeros(2,count);
    constreint = zeros(1,params.n+1);
    linsp = linspace (params.T, params.t0, count);
    solM=ode15s(@odeAM,linsp,[0,0],options,params);  
    m1=deval(solM, linsp);
    y=deval(soly, linsp);
    for j=params.n+1:-1:1
        params.k=j;
        index=count;
        for i=params.T:-(params.T-params.t0)/(count-1):params.t0
            dfdb([1, 2], index)=dudbi(i,params,y(:,index));
            index=index-1;
        end    
        constreint(params.n-j+2) = -trapz(linsp', dfdb(1,:).*m1(1,:)+dfdb(2,:).*m1(2,:));
    end
end



