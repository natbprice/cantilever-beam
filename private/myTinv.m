function U_norm=myTinv(U_stand_norm,opts)
% Transform from a standard gaussian space into design space
%
% Inputs:
% U_stand_norm  = vector of random variables in standard normal space
%
% Outputs:
% U_norm = vector of random variables in normalized design space
%==========================================================================

% Define distributions
name={'norm','norm'};
A=[500,1000];
B=[0.2,0.1].*A;

% Convert from standard normal space to design space
U=zeros(size(U_stand_norm));
U(:,1)=icdf(name{1},normcdf(U_stand_norm(:,1)),A(1),B(1));
U(:,2)=icdf(name{2},normcdf(U_stand_norm(:,2)),A(2),B(2));

% Convert from design space to normalized design space
U_norm=norm_design(U,'u',1,opts);

end

