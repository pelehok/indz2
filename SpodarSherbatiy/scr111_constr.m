clc
clear
p=[0.9,  -0.9,  -0.4,  -0.7,  0.5,  0.4,  0.001,  -0.1];
U0=p(6); %initial value of control function
yd=[7, 0.5]; %model parameters
t0=0; Tend=10; %time interval
y0=[5, 0.3]; %initial value
yy=[8, 0.6]; % y+
Umin=U0-0.1*abs(U0); Umax=U0+0.1*abs(U0);
z0=[0,0];

%Control function - piecewise-linear function
typeU='linear'
n=4; % - number of subdomains
n1=n+1; % - number of points
b0=U0*ones(1,n1) 

%Creation the structure of model parameters
params.p=p;
params.yd=yd;
params.t0=t0;params.Tend=Tend;
params.y0=y0;
params.yy=yy;
params.n=n;
params.U=b0;
params.lts=(Tend-t0)/n;
params.typeU=typeU;
params.j=6;

%Solving of direct problem for initial values of optimization parameters
options = odeset('RelTol',1e-9,'AbsTol',1e-9); 
sol=ode15s(@ode1_constr,[t0,Tend],y0,options,params);

params.sol=sol;
params.t=sol.x;
params.y=sol.y;
params.SAcoef = 0*ones(1,n1);

%==============================================================
% finding dPsi/db
%==============================================================


optcr = fun_optcr1(b0,params);
disp('Differences:');
for k=1:n1
bnew = b0;
step = 0.01*bnew(k);
bnew(k)=step+bnew(k);
optcr_i = fun_optcr1(bnew,params);
dpsidbi(k)=(optcr_i-optcr)/step;
end
dpsidbi


disp('DDM:');
dg1dy=Find_dg1dy (params.y);
SAcoef = 0*ones(1,n1);
for k=1:n1
    params.k=k;
    solz=ode15s(@odeZ,[t0, Tend],z0,options,params);
    linsp = linspace (t0, Tend, 100);
    zi=deval(solz, linsp);
    yy=deval(sol, linsp);
    dg1dy=Find_dg1dy(yy);
    SAcoef(k) = trapz(linsp', dg1dy(:,2).*zi(2,:)');
end
SAcoef
params.SAcoef = SAcoef;  

%================================================================
%Solution of optimization problem 
%================================================================

disp('Finding solution of optimization problem (Grad OFF)');
options1 = optimset('GradObj','off','MaxFunEvals',300);
un=Umin*ones(1,n1);
uv=Umax*ones(1,n1);
disp('Without Psi1:');
[x,FVAL] = fmincon(@fun_optcr1,b0,...
        [],[],[],[],un,uv,[],options1,params) % - Without Psi1 
% disp('All constraints:');
% [x,FVAL] = fmincon(@fun_optcr1,b0,...
%          [],[],[],[],un,uv,@constrPsi1,options1,params) % - all 

%=====================================================================

disp('Finding solution of optimization problem (Grad ON)');
options1 = optimset('GradObj','on','MaxFunEvals',300);
un=Umin*ones(1,n1);
uv=Umax*ones(1,n1);
disp('Without Psi1:');
[x,FVAL] = fmincon(@fun_optcr1,b0,...
        [],[],[],[],un,uv,[],options1,params) % - Without Psi1 
% disp('All constraints:');
% [x,FVAL] = fmincon(@fun_optcr1,b0,...
%          [],[],[],[],un,uv,@constrPsi1,options1,params) % - all 


%Calculation and plot of y for optimal value of optimization parameters
params.U=x;
[t,y]=ode15s(@ode1_constr,[t0,Tend],y0,options,params);
figure('Position',[100 100 1000 300])
subplot(1,2,1)
plot(t,y(:,1))
grid on
title('y1, optimal parameters,''fmincon''')
ylabel('y1')
xlabel('time, t')
subplot(1,2,2)
plot(t,y(:,2))
grid on
title('y2, optimal parameters,''fmincon''')
ylabel('y2')
xlabel('time, t')

%Calculation the value of optimization criteria and constraint
%for optimal value of optimization parameters
disp('Optim. criteria and constr. for optimal value of optim. parameters:');
[optcr] = fun_optcr1(x,params)
[c,constr] = constrPsi1(x,params)
disp('END finding solution of optimization problem');

%plot figures
tU=linspace(t0,Tend,300);
U_tU= Values_U_tU(tU,params);
figure
% p1=plot(tU,U_tU,'Color','red');
% hold on
params.U=b0;
U_tU= Values_U_tU(tU,params);
p2=plot(tU,U_tU,'Color','blue');
grid on
title('y, optimal parameters')
ylabel('y')
xlabel('time, t')
options1 = optimset('MaxFunEvals',1000);

%find optimal values of optimization parameters (without constraints)
disp('Optimal values of optimization parameters:')
[x,fval,exitflag,output] = fminsearch(@fun_optcr1,b0,options1,params)
params.U=x;
U_tU= Values_U_tU(tU,params);
hold on
p3=plot(tU,U_tU,'--','Color','black');
legend('Uinit','UoptSearch')

% solution of direct problem for optimal values of optimization parameters
[t,y]=ode15s(@ode1_constr,[t0,Tend],y0,options,params);
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

disp('Optimization criteria and constraint psi1 (without constraints):)')
[optcr] = fun_optcr1(x,params)
[c,constr] = constrPsi1(x,params)