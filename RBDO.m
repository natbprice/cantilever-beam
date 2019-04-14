function [xmin,fmin,counteval,stopflag,out,bestever]=RBDO(E,Opts,model_fidelity)
% Perform RBDO using CMA-ES
%
% Inputs:
% none
%
% Outputs:
% xmin = xmin from CMA-ES
% fmin = fmin from CMA-ES
% counteval = counteval from CMA-ES
% stopflag = stopflag from CMA-ES
% out = out from CMA-ES
% bestever = bestever from CMA-ES
%==========================================================================

P{1}=Opts;
P{2}=E;
P{3}=model_fidelity;

%% Bounds on design variable
opts.LBounds=zeros(2,1);        % Lower bounds on design variable
opts.UBounds=ones(2,1);         % Upper bounds on design variable

%% Convergence options
opts.MaxFunEvals=40e3;          % Maximum number of function evaluations
opts.TolFun=1e-2;               % Stop if function change is less than TolFun
opts.TolX=1e-5;                 % Stop if change in design variable is less than TolX

%% Display / data options
opts.DispFinal='on';            % Display messages like initial and final message
opts.DispModulo=1;              % Display messages after every i-th iteration
opts.LogPlot='on';              % Plot while running using output data files
opts.SaveFilename = ...         % Filename to save all variables
    'rbdo_variables_cmaes.mat';
opts.LogFilenamePrefix = ...    % Filename prefix for output data files
    'rbdo_outcmaes_';

%% Populations / restart options
% opts.PopSize=50;                % Initial population size
opts.Restarts=5;                % Number of restarts
opts.IncPopSize=2;              % Multiplier for population size before each restart

%% Noise options
opts.Noise.on=1;                % Turn on uncertainty handling
% opts.Noise.reevals=5;           % Number of re-evaluations for estimating uncertainty

%% Parallel options
opts.EvalParallel='yes';

%% Resume options
opts.Resume='no';               % Resume former run from saved variable file
if strcmp(opts.Resume,'yes')
    figure(324)
end

%% Call CMA-ES
start=rand(2,1);

% Run CMA-ES on MDA.m
[xmin,fmin,counteval,stopflag,out,bestever]=cmaes('RBDO_wrapper', start, 0.35, opts, P);


end