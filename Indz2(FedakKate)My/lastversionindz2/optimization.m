function b = optimization(params,s,yConstraint,bConstraint)
    if strcmp(s,'FDM')
      opt=@optFDM;
    end
    if strcmp(s,'DDM')
      opt=@optDDM; 
    end
    if strcmp(s,'AM')
      opt=@optAM; 
    end
    nonlinearconstraint=[];
    if yConstraint
        nonlinearconstraint=@constraintopt;
    end
    options = optimset('GradObj','on','MaxFunEvals',300);
    if bConstraint
        bLower=params.bLower*ones(1,params.n+1);
        bUpper=params.bUpper*ones(1,params.n+1);
    end
    if yConstraint||bConstraint
        [xFDM,FVAL,exitflag,output,lambda,grad] = fmincon(opt,params.b0,...
        [],[],[],[],bLower,bUpper,nonlinearconstraint,options,params);
    else
       [xFDM,FVAL,exitflag,output,lambda,grad] = fmincon(opt,params.b0,...
        [],[],[],[],bLower,bUpper,[],options,params);
    end
    b=xFDM
    FVAL
    grad
end

