function [E]=cond_sim_DoE(z,E,opts)
% Generate conditional simulations over a DoE
%
% Inputs:
% z = DoE in joint design / aleatory space
% m = number of conditional simulations
% opts = structure with inputs and options
%
% Outputs:
% E = updated error model structure with conditional simulation
%==========================================================================

%% Unpack variables
model=E.model;                              % STK error model
Zi=E.Zi{1};                                 % Initial DoE
Ei=E.Ei{1};                                 % Evaluations at initial DoE

%% Setting i=1 corresponds to DoE of conditional simulation
i=1;
    
%% Generate 1 conditional simulation over DoE z
e_rand=stk_generate_samplepaths(model,Zi,Ei,z,1);

%% Update error model structure with new DoE
kreq=stk_kreq_qr(model,[Zi;z]);
E.Zi{i+1}=[Zi;z];
E.Ei{i+1}=[Ei;e_rand(:,i)];
E.kreq{i+1}=kreq;

end

