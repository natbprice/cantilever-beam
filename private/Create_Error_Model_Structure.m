% Define structure with STK parameters for Kriging error model

E=struct();

%% Load data files
load('model')                               % Surrogate model of error in limit-state function

%% Define model structure
E.model=model;                              % STK model of model error
E.Zi=repmat({xi},2,1);                      % Initial DoE for model
E.Ei=repmat({zi},2,1);                      % Evaluations of true error on DoE
kreq = stk_kreq_qr(model, xi);
E.kreq=repmat({kreq},2,1);                  % Precomputed kreq object for speed

%% Initialize options that will be set internally
E.returnCalib=false;                        % Return calibrated model (set internally)
E.i=0;                                      % Index different DoE's (set internally)
E.returnKstd=false;                         % Return mean plus k standard devation offset
