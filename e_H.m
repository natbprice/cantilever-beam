function [e_out]=e_H(x,u,E,opts)
% Evaluate conditional simulation of Kriging error model or calibrated model
%
% Inputs:
% x = vector of design variables
% u = vector of aleatory variables
% E = structure with error model
% opts = structure with inputs and options
%
% Outputs:
% e_out = error prediction function evaluated at inputs
%==========================================================================

%% Unpack variables
returnCalib=E.returnCalib;                  % If true, return calibrated mean prediction
model=E.model;                              % STK error model
i=E.i;                                      % Current epistemic realization
Zi=E.Zi{i+1};                               % DoE
Ei=E.Ei{i+1};                               % Evaluations at DoE
kreq=E.kreq{i+1};                           % STK structure for fast evaluations
p_ini=opts.error.p_ini;                     % Number of samples in initial DoE

%% Return calibrated model or conditional simulation
z=[x,u];
if returnCalib
    e=stk_predict(model,Zi(1:p_ini+1,:),Ei(1:p_ini+1,:),z);
else
    e=stk_predict(model,{Zi,kreq},Ei,z);
end
if E.returnKstd
    e_out=e.mean+E.k*sqrt(e.var);
else
    e_out=e.mean;
end
    
end

