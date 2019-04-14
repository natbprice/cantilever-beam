% PRIVATE
%
% The files in this folder are specific to the cantilever beam problem
%
% Files
%   f                            - Objective function
%   g_H                          - Prediction of the limit-state function including Kriging model error
%   g_H_of_E                     - Prediction of the limit-state function as function of error realization
%   myTinv                       - Transform from a standard gaussian space into design space
%   norm_design                  - Transform from normalized space to design space (or vice versa)
%   Create_Error_Model_Structure - Define structure with STK parameters for Kriging error model
%   Create_Options_Structure     - Define structure will all options / parameters
%   g_H_true                     - High-fidelity limit-state function (Timoshenko model)
%   eul_def                      - Euler-Bernoulli model of tip deflection for cantilever beam
%   tim_def                      - Timoshenko model of tip deflection for cantilever beam
