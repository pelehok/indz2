function [Psi0,Psi1] = AM(params)
count = 100;
options = odeset('RelTol',1e-7,'AbsTol',1e-7);
dfdb=zeros(2,count);
SA=zeros(1,params.n+1);
constreint = zeros(1,params.n+1);
for j=params.n+1:-1:1
    params.k=j;
    linsp = linspace (params.T,params.t0, count);
    solAM0=ode15s(@odeAM,linsp,[0 0],options,params);
    solAM1=ode15s(@odeAM1,linsp,[0 0],options,params);
    m0=deval(solAM0, linsp);
    m1=deval(solAM1, linsp);
    y=deval(params.sol, linsp);
    index=count;
    for i=params.T:-(params.T-params.t0)/(count-1):params.t0
    dfdb([1, 2], index)=dudbi(i,params,y(:,index));
    index=index-1;
    end    
    SA(params.n+1-j+1) = -trapz(linsp', dfdb(1,:).*m0(1,:)+dfdb(2,:).*m0(2,:));
    constreint(params.n-j+2) = -trapz(linsp', dfdb(1,:).*m1(1,:)+dfdb(2,:).*m1(2,:));
end
Psi0=SA;
Psi1=constreint;
end

