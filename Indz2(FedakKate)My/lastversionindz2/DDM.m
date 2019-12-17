function [Psi0,Psi1] = DDM( params )
options = odeset('RelTol',1e-7,'AbsTol',1e-7); 
SA = zeros(1,params.n+1);
constreint = zeros(1,params.n+1);
for j=1:params.n+1
    params.k=j;
    solDDM=ode15s(@odeDDM,[params.t0,params.T],[0 0],options,params);
    linsp = linspace (params.t0, params.T, 100);
    zi=deval(solDDM, linsp);
    y=deval(params.sol, linsp);
    g0=[(zeros(1,size(y,2)))', 2*y(2,:)'-1];
    g1=[2*(abs(y(1,:)-params.yMax(1))+y(1,:)-params.yMax(1))'.*...
    (sign(y(1,:)-params.yMax(1))+1)', (zeros(1,size(y,2)))'];
    constreint(j) = trapz(linsp',g1(:,1).*zi(1,:)');
    SA(j) = trapz(linsp', g0(:,2).*zi(2,:)');
end
Psi0=SA;
Psi1=constreint;
end

