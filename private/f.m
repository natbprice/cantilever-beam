function [m] = f(x_norm,u_norm,E,opts)
% Objective function
%
% Inputs:
% x  = vector of design variables
% u = vector of aleatory variables
% E = structure with error model
% opts = structure with inputs and options
%
% Outputs:
% ff = cost function evaluated as inputs
%==========================================================================
% Scale normalized variable
x=norm_design(x_norm,'x',0,opts);

b=x(:,1);           % Cross section base width (in)
h=x(:,2);           % Cross section height (in)

m=b.*h;

end