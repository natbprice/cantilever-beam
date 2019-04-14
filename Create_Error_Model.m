function [model,xi,zi]=Create_Error_Model(opts)
% Create Kriging model for error between beam models
%
% Inputs:
% opts = structure with inputs and options
%
% Outputs:
% model = STK Kriging model
% xi = DoE
% zi = Evaluations on DoE
%==========================================================================

%% Design of experiment
% Corner points
lb=ones(1,4);
ub=[2,2,2,2];
xff=fullfact([2,2,2,2]);
xff_norm=bsxfun(@rdivide,bsxfun(@minus,xff,lb),(ub-lb));
xi=xff_norm;

%% Evaluate model
e=@(x,u) g_H_true(x,u,opts)-g_H_of_E(x,u,0,opts);
zi=e(xi(:,1:2),xi(:,3:4));

%% Fit GP model
% Set covariance function and model order
dim=4;
model=stk_model('stk_gausscov_aniso',dim);

% Compute an initial guess for the covariance parameters and a reasonable
% log-variance for a small "regularization noise"
box = [zeros(1,dim);ones(1,dim)];
[param0, model.lognoisevariance]=stk_param_init(model,xi,zi,box);

% Estimate covariance parameters from data
model.param=stk_param_estim(model,xi,zi,param0);
model.lognoisevariance=log(eps);

end