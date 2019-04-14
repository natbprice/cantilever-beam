function [xopt,fopt] = Design_Opt_MCS(n,E,opts)
% MCS of deterministic design optimization for different error realizations
%
% Inputs:
% n = number of MCS samples
% E = structure with error model
% opts = structure with inputs and options
%
% Outputs:
% xopt  = matrix of optimum design variables [n x dim]
% fopt = matrix of objective function [n x dim]
%==========================================================================

udet=opts.design.udet;
DoE=opts.error.DoE_for_cond_sims;

Ein=E;
Ein.i=1;
xopt=zeros(n,5);
fopt=zeros(n,1);
parfor i=1:n
    E=cond_sim_DoE(DoE,Ein,opts);
    E.returnCalib=false;
    [xopt(i,:),fopt(i,:)]=Design_Opt(0,udet,E,opts,rand(1,5));
end

end