function dmiu = odeMiu(t,miu,params)

% yii=interp1(params.t, params.y, t);
yi=deval (params.sol, t);

dfdy=Grads(params, yi);
dg1dy = Find_dg1dy(yi);

dmiu=-dfdy' * miu - dg1dy';

end
