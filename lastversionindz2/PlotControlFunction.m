function [] = PlotControlFunction(params,x)
tU=linspace(params.t0,params.T,300);
U_tU= du(tU,params);
figure
params.U=params.b0;
U_tU= du(tU,params);
p2=plot(tU,U_tU,'Color','blue');
grid on
title('U, control function')
ylabel('U')
xlabel('time, t')
params.U=x;
U_tU= du(tU,params);
hold on
p3=plot(tU,U_tU,'--','Color','black');
legend('Uinit','UoptSearch')
end

