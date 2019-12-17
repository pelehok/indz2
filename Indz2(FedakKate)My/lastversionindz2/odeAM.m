function dm = odeAM( t,m,params )
y=deval (params.sol,t);
df1dy1= params.p(1) - (params.p(4)*y(2)*(exp(-params.p(5)*y(1)) - 1))/(params.p(6) + params.p(7)*y(2)) +...
    (params.p(4)*params.p(5)*y(1)*y(2)*exp(-params.p(5)*y(1)))/(params.p(6) + params.p(7)*y(2));
df1dy2= (params.p(4)*params.p(7)*y(1)*y(2)*(exp(-params.p(5)*y(1)) - 1))/(params.p(6) + params.p(7)*y(2))^2 -...
    (params.p(4)*y(1)*(exp(-params.p(5)*y(1)) - 1))/(params.p(6) + params.p(7)*y(2));
df2dy1= (params.p(4)*params.p(5)*params.p(8)*y(1)*y(2)*exp(-params.p(5)*y(1)))/(params.p(6) + params.p(7)*y(2)) -...
    (params.p(4)*params.p(8)*y(2)*(exp(-params.p(5)*y(1)) - 1))/(params.p(6) + params.p(7)*y(2));
df2dy2= params.p(2) + 2*params.p(3)*y(2) - (params.p(4)*params.p(8)*y(1)*(exp(-params.p(5)*y(1)) - 1))/...
    (params.p(6) + params.p(7)*y(2)) + (params.p(4)*params.p(7)*params.p(8)*y(1)*y(2)*(exp(-params.p(5)*y(1)) - 1))/...
    (params.p(6) + params.p(7)*y(2))^2;
dfdy=[df1dy1, df1dy2;
      df2dy1, df2dy2];
g0=[(zeros(1,size(y,2)))', 2*y(1,:)'-1];
dm=-dfdy' * m - g0';
end

