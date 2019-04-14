function [xopt,fopt] = Design_Opt(k,udet,E,opts,x0)
% Perform deterministic safety-margin-based design optimization
%
% Inputs:
% k = standard deviation offset
% udet = vector of conservative values used in place of aleatory variables
% E = structure with error model
% opts = structure with inputs and options
% x0 = starting point for optimization (optional)
%
% Outputs:
% xopt  = vector of optimum design variables
% fopt = cost function evaluated at optimum
%==========================================================================

%% Specify dimensions of optimization problem
dim=opts.design.dim;

%% Define default starting point
if nargin<5
    x0=0.5*ones(1,dim);
end

%% Define objective function and constraints
% Objective function
function F=myObjFun(x_norm)
    F=f(x_norm,udet,E,opts);
end

% Constraint function
function [c,ceq]=myConFun(x_norm)
    E.returnKstd=true;
    E.k=-k;
    g1=g_H(x_norm,udet,E,opts);
    c=-g1;
    ceq=[];
end

%% Call fmincon to perform design optimization
% Define options for fmincon
options = optimoptions('fmincon','Algorithm','sqp','Display',opts.design.display,'TolCon',1e-4,'TolFun',1e-4,'TolX',1e-4);

% Call fmincon
[xopt,fopt,exitflag,output]=fmincon(@myObjFun,x0,[],[],[],[],zeros(1,dim),ones(1,dim),@myConFun,options);

% Try random restarts if optimization fails
num_restarts=10;
if exitflag<=0
    j=1;
    while exitflag<=0 && j<=num_restarts
        [xopt,fopt,exitflag,output]=fmincon(@myObjFun,rand(1,dim),[],[],[],[],zeros(1,dim),ones(1,dim),@myConFun,options);
        j=j+1;
    end
    if exitflag<=0
        display('Design optimization failed')
    end
end

end


