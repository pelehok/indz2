function sol = DirectProblem(params)
options = odeset('RelTol',1e-7,'AbsTol',1e-7); 
sol=ode15s(@ode,[params.t0,params.T],params.y0,options,params);
end

