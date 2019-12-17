clc;
clear;
p=[0.9,  -0.9,  -0.4,  -0.7,  0.5,  0.4,  0.001,  -0.1];
U0=p(8); 
yd=[7, 0.5]; 
t0=0;
T=10;
y0=[5, 0.3]; 
yMax=[8, 0.6];
bLower=U0-0.1*abs(U0);
bUpper=U0+0.1*abs(U0);
n=1; 
npoint=n+1;
b0=U0*ones(1,npoint);
params.p=p;
params.yd=yd;
params.bLower=bLower;
params.bUpper=bUpper;
params.t0=t0;
params.T=T;
params.y0=y0;
params.yMax=yMax;
params.n=n;
params.U=b0;
params.b0=b0;
params.step=(T-t0)/n;
params.j=8;
%Slove direct problem
sol=DirectProblem(params);
params.sol=sol;
params.t=sol.x;
params.y=sol.y;
params.SA = zeros(1,npoint);
%Plot(params);
%Find gradien throught fmincon with 0 iteration
tic;
GradFmincon(npoint,b0,params,bLower,bUpper);
toc;
%FDM 
disp('FDM:');
tic;
[Psi0_fdm, Psi1_fdm] = FDM(b0,params)
toc;
%DDM
disp('DDM:');
tic;
[Psi0_ddm, Psi1_ddm]=DDM(params)
toc;
% %AM
disp('AM:');
tic;
[Psi0_am, Psi1_am]=AM(params)
toc;
%----------------------------------------------
%Optimization 
% tic;
% bnew=optimization(params,'FDM',true,true) ;
% toc;  
% tic;
%  bnew2=optimization(params,'DDM',true,true) ;
%  toc;
% tic;
% bnew3=optimization(params,'AM',true,true) ;
% toc;
%-----------------------------------------------
%disp('PlotOptimal');
%PlotOptimal(params,bnew)
%-----------------------------------------------
%find optimal values of optimization parameters (without constraints)
%disp('Optimal values of optimization parameters:')
%bopt=OptimalOptimization(params);
%-------------------------------------------------------
% plot control function 
%PlotControlFunction(params,bopt);