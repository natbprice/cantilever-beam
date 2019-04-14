% Define structure will all options / parameters

opts=struct();

%% Load data files
load('srgt_models')                         % Surrogate model of expected objective function as function of safety margins
load('DoE_for_cond_sims')                   % DoE for conditional simulations
load('model')                               % Surrogate model of error in limit-state function

%% Safety margin optimization
opts.margins.method='safety';               % Redesign method
opts.margins.numRestarts=5;                 % Number of restarts for CMA-ES
opts.margins.lb=[0,1e-3,1e-3,0];            % Lower bound on safety margins
opts.margins.ub=[4,4,4,4];                  % Upper bound on safety margins
opts.margins.alpha_bar=0.05;                % Target confidence level for reliability constraint
opts.margins.pf_bar=normcdf(-3);            % Target mean probability of failure
opts.margins.pre_bar=0.50;                  % Target probability of redesign
opts.margins.w1=1e3;                        % Weight on probability of failure constraint
opts.margins.w2=1e3;                        % Weight on probability of redesign constraint
opts.margins.SaveFilename=...               % Filename for CMA-ES variables
    'test_vars.mat';
opts.margins.LogFilenamePrefix=...          % Filename prefix for CMA-ES logs
    'test_';
opts.margins.Srgts=srgts;                   % Surrogates for fast safety margion optimization

%% Epistemic uncertainty definition
opts.error.display=false;                    % Display results for MCS
opts.error.mstart=500;                      % Default number of MCS samples
opts.error.mmax=1e4;                        % Maximum number of MCS samples
opts.error.cov_target=[5,0.5];             % Target COV's (%) for alpha and Ef                 
opts.error.DoE_for_cond_sims=n_norm;        % Define DoE for performing conditional simulations
opts.error.p_ini=length(zi);                % Define number of initial samples in DoE

%% Figure options
opts.fig.createHists=false;                  % Create histograms (MCS only)
opts.fig.saveFiles=true;                    % Export figures to file
opts.fig.fname='_perf';                     % Suffix of image filenames
opts.fig.plot_redesign=false;               % Plot true redesign outcome
opts.fig.xlabels={...                       % Design variable labels
      'Thickness (in)','Height (in)'};
opts.fig.flabel={...                        % Objective function label
    'Area of cross section (in^2)'};

% Bounds on axes of figures (optional)
opts.fig.xlims.f=[9,10.6];
opts.fig.xlims.x=[3.3,3.65;2.65,2.9];
opts.fig.xlims.margins=[-4e-4,6e-4];
opts.fig.xlims.beta=[1,6];
opts.fig.xlims.mpp=[600,1000;1100,1300];

opts.fig.ylims.f=[0,2e3];
opts.fig.ylims.margins=[0,600];
opts.fig.ylims.beta=[0,600];
opts.fig.ylims.mpp=[0,110];

opts.fig.caxis.mpp=[0,100];

%% Reliability analysis
opts.reliability.calcPf=true;               % Option to not calculate probability of failure
opts.reliability.method='form';             % Method for calculating probability of failure
opts.reliability.vectorial=true;            % Specify if limit-state is vecorized
opts.reliability.dim=2;                     % Number of aleatory variables
opts.reliability.lb=[300,800];              % Lower bound on aleatory variables
opts.reliability.ub=[1200,1700];            % Upper bound on aleatory variables

% FORM options
opts.reliability.form.eps=1e-4;             % Tolerance for optimization
opts.reliability.form.solver='sqp';         % Optimization algorithm
opts.reliability.form.display='none';       % Display option for CODES toolbox ('none', 'final', 'iter')

% MCS options
opts.reliability.mcs.nstart=1e6;            % Starting sample size for MCS
opts.reliability.mcs.nmax=1e10;             % Maximum sample size for MCS
opts.reliability.mcs.cov=5;                 % Target COV (%) for MCS estimate
opts.reliability.mcs.display=false;         % Display option for CODES toolbox (true/false)

%% Deterministic design optimization
opts.design.dim=2;                          % Number of design variables
opts.design.display='none';                 % Display option for fmincon
opts.design.numRestarts=10;                 % Number of restarts if fmincon fails
opts.design.lb=[2.5,1.5];                   % Lower bound on design variables
opts.design.ub=[5.5,4.5];                   % Upper bound on design variables
opts.design.udet=[0.494132779173301,...     % Conservative values for aleatory variables
    0.415032388338957];

%% Objective function MCS
opts.obj.nstart=50;                         % Starting sample size for MCS
opts.obj.nmax=1e3;                          % Maximum sample size for MCS
opts.obj.vectorial=false;                   % Limit-state function is vectorized for MCS
opts.obj.cov=0.5;                           % Target COV (%) for MCS estimate

