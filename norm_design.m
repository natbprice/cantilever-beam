function [X_out] = norm_design(X_in,x_or_u,flag,opts)
% Transform from normalized space to design space (or vice versa)
%
% Inputs:
% X_in  = input design variables
% opts = structure with inputs and options
%
% Outputs:
% X_out = output design variables
%==========================================================================

% Unpack options
dim_x=opts.design.dim;
dim_u=opts.reliability.dim;

% Bounds on design variables
lb=[opts.design.lb,opts.reliability.lb];
ub=[opts.design.ub,opts.reliability.ub];

% Input is in x-space (design) or u-space (aleatory)
switch lower(x_or_u)
    case 'x'
        lb=lb(1:dim_x);
        ub=ub(1:dim_x);
    case 'u'
        lb=lb(dim_x+1:dim_x+dim_u);
        ub=ub(dim_x+1:dim_x+dim_u);
    otherwise
        lb=lb(str2double(x_or_u));
        ub=ub(str2double(x_or_u));    
end

% Convert normalized to design space or design space to normalized
if flag==0
    % Normalized space to design space
    X_out=bsxfun(@plus,lb,bsxfun(@times,X_in,(ub-lb)));
else  
    % Design space to normalized space
    X_out=bsxfun(@rdivide,bsxfun(@minus,X_in,lb),(ub-lb));
end

end

