function dz = odeZ(t,z,params)

% yii=interp1(params.t, params.y, t);
yi=deval (params.sol, t);

dfdb=find_dfdb(t,params,yi);
dfdy=Grads(params, yi);

dz=dfdb + dfdy * z;

end