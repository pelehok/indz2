clc
clear
p=[0.9,  -0.9,  -0.4,  -0.7,  0.5,  0.4,  0.001,  -0.1];
U0=p(6); %initial value of control function
yd=[7, 0.5]; %model parameters
t0=0; Tend=10; %time interval
y0=[5, 0.3]; %initial value
yy=[8, 0.6]; % y+

%Solving of direct problem for initial values of optimization parameters
options = odeset('RelTol',1e-4,'AbsTol',1e-5); 
[t,y]=ode15s(@ode1,[t0,Tend],y0,options,p);

%plot graphs
figure('Position',[100 100 1000 300])
subplot(1,2,1)
plot(t,y(:,1))
grid on
title('y1, initial parameters')
ylabel('y1')
xlabel('time, t')
subplot(1,2,2)
plot(t,y(:,2))
grid on
title('y2, initial parameters')
ylabel('y2')
xlabel('time, t')
figure
plot(y(:,1),y(:,2))
grid on
title('Phase portrait')
ylabel('y2')
xlabel('y1')