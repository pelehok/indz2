function [ dPsi0dbi, dPsi1dbi ] = Differences(u, params )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
params.U=u;
optcr = fun_optcr1(params.U,params);
[c,constr] = constrPsi1(params.U, params);
for k=1:params.n+1
bnew = params.U;
step = 0.00001*bnew(k);
bnew(k)=step+bnew(k);
optcr_i = fun_optcr1(bnew,params);
[c,constr_i]=constrPsi1(bnew,params);
dPsi0dbi(k)=(optcr_i-optcr)/step;
dPsi1dbi(k)=(constr_i-constr)/step;
end

end

