function [dPsi0dbi, dPsi1dbi] = FDM(u, params )
params.U=u;
criter = criteria(params.U,params);
constr = constraint(params.U, params);
for k=1:params.n+1
bnew = params.U;
step = 0.00001*bnew(k);
bnew(k)=step+bnew(k);
criter_i = criteria(bnew,params);
constr_i=constraint(bnew,params);
dPsi0dbi(k)=(criter_i-criter)/step;
dPsi1dbi(k)=(constr_i-constr)/step;
end

end

