function [gPsi0,gPsi1] =GradFmincon( npoint,b0,params,Umin,Umax)
options1 = optimset('GradObj','off','MaxFunEvals',0);
un=Umin*ones(1,npoint);
uv=Umax*ones(1,npoint);
[x,fval,exitflag,output,lambda,Psi0] = fmincon(@criteria,b0,[],[],[],[],...
    un,uv,[],options1,params);
[x,fval,exitflag,output,lambda,Psi1] = fmincon(@constraint,b0,[],[],[],[],...
    un,uv,[],options1,params);
gPsi0=Psi0
gPsi1=Psi1
end

