clc
clear
p=[0.9,  -0.9,  -0.4,  -0.7,  0.5,  0.4,  0.001,  -0.1];
j=4;
U0=p(j); %initial value of control function
yd=[7, 0.5]; %model parameters
t0=0; Tend=10; %time interval
y0=[5, 0.3]; %initial value
yy=[8, 0.6]; % y+
Umin=U0-0.11*abs(U0); Umax=U0+0.11*abs(U0);
z0=[0,0];
miu0=[0,0];

%Control function - piecewise-linear function
typeU='linear'
n=1; % - number of subdomains
n1=n+1; % - number of points
b0=U0*ones(1,n1) 

%Creation the structure of model parameters
params.p=p;
params.yd=yd;
params.t0=t0;params.Tend=Tend;
params.y0=y0; params.z0=z0; params.miu0=miu0;
params.yy=yy;
params.n=n;
params.U=b0; params.b0=b0;
params.lts=(Tend-t0)/n;
params.typeU=typeU;
params.j=j;

%Solving of direct problem for initial values of optimization parameters
options = odeset('RelTol',1e-7,'AbsTol',1e-7); 
sol=ode15s(@ode1_constr,[t0,Tend],y0,options,params);

params.sol=sol;
params.t=sol.x;
params.y=sol.y;
params.SAcoef = 0*ones(1,n1);

%==============================================================
% finding dPsi/db
%==============================================================

disp('Gradient by fmincon:')
options1 = optimset('GradObj','off','MaxFunEvals',0);
un=Umin*ones(1,n1);
uv=Umax*ones(1,n1);
[x,fval,exitflag,output,lambda,gradPsi0] = fmincon(@fun_optcr1,b0,[],[],[],[],...
    un,uv,[],options1,params);
[x,fval,exitflag,output,lambda,gradPsi1] = fmincon(@constrPsi111,b0,[],[],[],[],...
    un,uv,[],options1,params);
gradPsi0
gradPsi1
    
%===============================================

disp('Differences:');
[dPsi0dbi, dPsi1dbi] = Differences(b0,params)

%================================================

disp('DDM:');
dg1dy=Find_dg1dy (params.y);
SAcoefDDM = 0*ones(1,n1);
for k=1:n1
    params.k=k;
    solz=ode15s(@odeZ,[t0, Tend],z0,options,params);
    linsp = linspace (t0, Tend, 100);
    zi=deval(solz, linsp);
    yy=deval(sol, linsp);
    dg1dy=Find_dg1dy(yy);
    SAcoefDDM(k) = trapz(linsp', dg1dy(:,2).*zi(2,:)');
end
SAcoefDDM
params.SAcoefDDM = SAcoefDDM;  

[c,constr, gradPsi1DDM] = constrPsi1DDM(params.U, params);
gradPsi1DDM

%=================================================

disp('AM:');
count = 100;
dfdb=0*ones(2,count);
SAcoefAM = 0*ones(1,n1);
for k=n1:-1:1
    params.k=k;
    linsp = linspace (Tend, t0, count);
    solmiu=ode15s(@odeMiu,linsp,miu0,options,params);  
    miu=deval(solmiu, linsp);
    yy=deval(sol, linsp);
    index=count;
    for q=Tend:-(Tend-t0)/(count-1):t0
        dfdb([1, 2], index)=find_dfdb(q,params,yy(:,index));
        index=index-1;
    end    
    SAcoefAM(n1-k+1) = -trapz(linsp', dfdb(1,:).*miu(1,:)+dfdb(2,:).*miu(2,:));
end
SAcoefAM
params.SAcoefAM = SAcoefAM;  

[c,constr, gradPsi1AM] = constrPsi1AM(params.U, params);
gradPsi1AM

% %================================================================
% %Solution of optimization problem 
% %================================================================
% 
% disp('Finding solution of optimization problem (Grad ON, FDM)');
% options1 = optimset('GradObj','on','MaxFunEvals',300);
% un=Umin*ones(1,n1);
% uv=Umax*ones(1,n1);
% % disp('Without Psi1:');
% % t=cputime;
% % [xFDM,FVAL,exitflag,output,lambda,grad] = fmincon(@fun_optcrFDM,b0,...
% %         [],[],[],[],un,uv,[],options1,params) % - Without Psi1 
% % timeFDM = cputime-t    
% disp('All constraints:');
% t=cputime;
% [xFDM,FVAL,exitflag,output,lambda,grad] = fmincon(@fun_optcrFDM,b0,...
%          [],[],[],[],un,uv,@constrPsi1,options1,params) % - all 
% timeFDM = cputime-t  
% 
% %=====================================================================
% 
% disp('Finding solution of optimization problem (Grad ON, DDM)');
% options1 = optimset('GradObj','on','MaxFunEvals',300);
% un=Umin*ones(1,n1);
% uv=Umax*ones(1,n1);
% % disp('Without Psi1:');
% % t=cputime;
% % [xDDM,FVAL,exitflag,output,lambda,grad] = fmincon(@fun_optcrDDM,b0,...
% %         [],[],[],[],un,uv,[],options1,params) % - Without Psi1
% % timeDDM = cputime-t    
% disp('All constraints:');
% t=cputime;
% [xDDM,FVAL,exitflag,output,lambda,grad] = fmincon(@fun_optcrDDM,b0,...
%          [],[],[],[],un,uv,@constrPsi1,options1,params) % - all 
% timeDDM = cputime-t
% 
% %==============================================================
% 
% disp('Finding solution of optimization problem (Grad ON, AM)');
% options1 = optimset('GradObj','on','MaxFunEvals',300);
% un=Umin*ones(1,n1);
% uv=Umax*ones(1,n1);
% % disp('Without Psi1:');
% % t=cputime;
% % [xAM,FVAL,exitflag,output,lambda,grad] = fmincon(@fun_optcrAM,b0,...
% %         [],[],[],[],un,uv,[],options1,params) % - Without Psi1
% % timeAM = cputime-t    
% disp('All constraints:');
% t=cputime;
% [xAM,FVAL,exitflag,output,lambda,grad] = fmincon(@fun_optcrAM,b0,...
%          [],[],[],[],un,uv,@constrPsi1,options1,params) % - all 
% timeAM = cputime-t
% 
% %==================================================================
% %Calculation and plot of y for optimal value of optimization parameters
% %==================================================================
% 
% % params.U=xFDM;
% % [t,y]=ode15s(@ode1_constr,[t0,Tend],y0,options,params);
% % figure('Position',[100 100 1000 300])
% % subplot(1,2,1)
% % plot(t,y(:,1))
% % grid on
% % title('y1, optimal parameters,''fmincon FDM''')
% % ylabel('y1')
% % xlabel('time, t')
% % subplot(1,2,2)
% % plot(t,y(:,2))
% % grid on
% % title('y2, optimal parameters,''fmincon FDM''')
% % ylabel('y2')
% % xlabel('time, t')
% % 
% % params.U=xDDM;
% % [t,y]=ode15s(@ode1_constr,[t0,Tend],y0,options,params);
% % figure('Position',[100 100 1000 300])
% % subplot(1,2,1)
% % plot(t,y(:,1))
% % grid on
% % title('y1, optimal parameters,''fmincon DDM''')
% % ylabel('y1')
% % xlabel('time, t')
% % subplot(1,2,2)
% % plot(t,y(:,2))
% % grid on
% % title('y2, optimal parameters,''fmincon DDM''')
% % ylabel('y2')
% % xlabel('time, t')
% % 
% % params.U=xAM;
% % [t,y]=ode15s(@ode1_constr,[t0,Tend],y0,options,params);
% % figure('Position',[100 100 1000 300])
% % subplot(1,2,1)
% % plot(t,y(:,1))
% % grid on
% % title('y1, optimal parameters,''fmincon AM''')
% % ylabel('y1')
% % xlabel('time, t')
% % subplot(1,2,2)
% % plot(t,y(:,2))
% % grid on
% % title('y2, optimal parameters,''fmincon AM''')
% % ylabel('y2')
% % xlabel('time, t')
% 
% %===========================================================
% %find optimal values of optimization parameters (without constraints)
% %=============================================================
% 
% disp('Optimal values of optimization parameters:')
% options1 = optimset('MaxFunEvals',1000);
% [x,fval] = fminsearch(@fun_optcr1,b0,options1,params)
% params.U=x;
% 
% % solution of direct problem for optimal values of optimization parameters
% [t,y]=ode15s(@ode1_constr,[t0,Tend],y0,options,params);
% figure('Position',[100 100 1000 300])
% subplot(1,2,1)
% plot(t,y(:,1))
% grid on
% title('y1, optimal parameters,''fminsearch''')
% ylabel('y1')
% xlabel('time, t')
% subplot(1,2,2)
% plot(t,y(:,2))
% grid on
% title('y2, optimal parameters,''fminsearch''')
% ylabel('y2')
% xlabel('time, t')
% 
% %plot control function 
% tU=linspace(t0,Tend,300);
% U_tU= Values_U_tU(tU,params);
% figure
% params.U=b0;
% U_tU= Values_U_tU(tU,params);
% p2=plot(tU,U_tU,'Color','blue');
% grid on
% title('U, control function')
% ylabel('U')
% xlabel('time, t')
% params.U=x;
% U_tU= Values_U_tU(tU,params);
% hold on
% p3=plot(tU,U_tU,'--','Color','black');
% legend('Uinit','UoptSearch')
