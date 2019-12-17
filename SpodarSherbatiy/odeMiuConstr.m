function dmiu = odeMiuConstr(t,miu,params)

% yii=interp1(params.t, params.y, t);
yi=deval (params.sol, t);

dfdy=Grads(params, yi);
dg1dy = [2*(abs(yi(1,:)-params.yy(1))+yi(1,:)-params.yy(1))'.*...
    (sign(yi(1,:)-params.yy(1))+1)', (0*ones(1,size(yi,2)))'];

dmiu=-dfdy' * miu - dg1dy';

end
