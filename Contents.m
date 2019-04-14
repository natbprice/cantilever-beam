% CANTILVER_BEAM_EXAMPLE
%
% Cantilever beam bending example with Kriging error model and quantile
% reliability constraint.
%
% Private folder contains files specific to beam example.
%
% The Matlab code requires:
% (1) Computational Optimal Design of Engineering Systems (CODES) toolbox
% (2) Covariance Matrix Adaptation Evolution Strategy (CMA-ES) code
% (3) Matlab Statistics Toolbox
% (4) Small Toolbox for Kriging (STK)
% (5*) Matlab Parallel Computing Toolbox is recommended
%
% Files
%   Design_Opt                 - Perform deterministic safety-margin-based design optimization
%   Reliability_Analysis       - Perform reliability analysis using CODES Toolbox
%   Safety_Margin_Opt          - Optimize safety margins (CMA-ES)
%   All_Tradeoff_Curves        - Create multiple tradeoff curves: safety, performance, mixed
%   Create_Tradeoff_Curve      - Create tradeoff curve for expected performance vs probability of redesign
%   Design_Process_Sim_Adapt   - Perform MCS of deterministic design process with adaptive sample size
%   Design_Process_Sim_srgt    - Call STK surrogate models in place of Design_Process_Sim_Adapt.m
%   e_H                        - Evaluate conditional simulation of Kriging error model or calibrated model
%   fig_export                 - Save figure to image file
%   norm_design                - Transform from normalized space to design space (or vice versa)
%   Objective_Function_MCS     - Calculate mean value of objective function w.r.t. aleatory uncertainty
%   VarRatio                   - Calculate ratio of epistemic to aleatory variance
%   Safety_Margin_DoE          - Create DoE for fitting surrogates as a function of safety margins
%   Create_Safety_Margin_Srgts - Create surrogates as a function of safety margins
%   cond_sim_DoE               - Generate conditional simulations over a DoE
%   Design_Opt_MCS             - MCS of deterministic design optimization for different error realizations
%   MCS_tradeoff_curves        - Perform MCS for each point on tradeoff curve
%   PlotHists                  - Generate histograms from MCS data
%   RBDO                       - Perform RBDO using CMA-ES
%   RBDO_wrapper               - Calculate penalized objective function considering probability of failure
%   Eval_xini                  - Evaluate true probability of failure and safety margin for tradeoff curve
%   Create_Error_Model         - Create Kriging model for error between beam models
%   plot_error                 - Create plots illustrating design optimization and reliability assesment
