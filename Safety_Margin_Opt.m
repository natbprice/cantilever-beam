function [xmin,fmin,counteval,stopflag,out,bestever] = Safety_Margin_Opt(P1)
% Optimize safety margins (CMA-ES)
%
% Inputs:
% P1 = structure with problem parameters / options
%
% Outputs:
% xmin = xmin from CMA-ES
% fmin = fmin from CMA-ES
% counteval = counteval from CMA-ES
% stopflag = stopflag from CMA-ES
% out = out from CMA-ES
% bestever = bestever from CMA-ES
%==========================================================================

%% Set number of design variables
if strcmp(P1.margins.method,'mixed')
    dim=4;
elseif strcmp(P1.margins.method,'safety') || strcmp(P1.margins.method,'performance')
    dim=3;
end

%% Bounds on design variable
opts.LBounds=zeros(dim,1);      % Lower bounds on design variable
opts.UBounds=ones(dim,1);       % Upper bounds on design variable

%% Convergence options
opts.MaxFunEvals=100e3;         % Maximum number of function evaluations
opts.TolFun=1e-6;               % Stop if function change is less than TolFun
opts.TolX=1e-6;                 % Stop if change in design variable is less than TolX

%% Display / data options
opts.DispFinal='on';            % Display messages like initial and final message
opts.DispModulo=1;              % Display messages after every i-th iteration
opts.LogPlot='on';              % Plot while running using output data files
opts.SaveFilename = ...         % Filename to save all variables
    P1.margins.SaveFilename;
opts.LogFilenamePrefix = ...    % Filename prefix for output data files
    P1.margins.LogFilenamePrefix;

%% Populations / restart options
% opts.PopSize=10*dim;                  % Initial population size
opts.Restarts=P1.margins.numRestarts;   % Number of restarts
opts.IncPopSize=2;                      % Multiplier for population size before each restart

%% Noise options
opts.Noise.on=0;                % Turn on uncertainty handling
% opts.Noise.reevals=3;         % Number of re-evaluations for estimating uncertainty

%% Parallel options
opts.EvalParallel='yes';        % Call objective function with multiple evaluations 

%% Resume options
opts.Resume='no';               % Resume former run from saved variable file

%% Call CMA-ES
start=rand(dim,1);
[xmin, fmin, counteval, stopflag, out, bestever]=cmaes('Design_Process_Sim_srgt', start, 0.35, opts, P1);

end

