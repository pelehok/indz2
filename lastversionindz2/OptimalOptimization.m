function [b] = OptimalOptimization(params)
options1 = optimset('MaxFunEvals',1000);
options = odeset('RelTol',1e-7,'AbsTol',1e-7);
[x,fval] = fminsearch(@criteria,params.b0,options1,params)
b=x;
params.U=x
[t,y]=ode15s(@ode,[params.t0,params.T],params.y0,options,params);
figure('Position',[100 100 1000 300])
subplot(1,2,1)
plot(t,y(:,1))
grid on
title('y1, optimal parameters,''fminsearch''')
ylabel('y1')
xlabel('time, t')
subplot(1,2,2)
plot(t,y(:,2))
grid on
title('y2, optimal parameters,''fminsearch''')
ylabel('y2')
xlabel('time, t')
end

