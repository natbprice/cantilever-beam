function [g]=g_H_of_E(x_norm,u_norm,e_i,opts)
% Prediction of the limit-state function as function of error realization
%
% Inputs:
% x = vector of design variables
% u = vector of aleatory variables
% e_i = realization of model error
% opts = structure with inputs and options
%
% Outputs:
% g_out = prediction of limit-state function evaluated at inputs
%==========================================================================

% Scale normalized variable
x=norm_design(x_norm,'x',0,opts);
u=norm_design(u_norm,'u',0,opts);

dtip=2.25*(0.1^3);          % Allowable tip displacement
b=x(:,1);                   % Cross section base width (in)
h=x(:,2);                   % Cross section height (in)
Px=u(:,1);                  % Tip load in x-direction (lbs)
Py=u(:,2);                  % Tip load in y-direction (lbs)
L=100*0.1;                  % Length (in)
E2=29e6;                    % Elastic modulus for steel (psi)
G=11.2e6;                   % Shear modulus for steel (psi)
x=0;                        % Distance from free end (in)
Iz=(b.*h.^3)./12;           % Area moment of inertia (in^4)
Iy=(h.*b.^3)./12;           % Area moment of inertia (in^4)

% Eq. xiii, pg. 59, Aircraft Structures for Engineering Students, 4th
% ed, Megson, 2007
dx=((Px*x^3)./(6*E2*Iz)) - ((Px*L^2*x)./(2*E2*Iz)) + ((Px*L^3)./(3*E2*Iz));
dy=((Py*x^3)./(6*E2*Iy)) - ((Py*L^2*x)./(2*E2*Iy)) + ((Py*L^3)./(3*E2*Iy));
d=sqrt(dx.^2+dy.^2);

g=dtip-d+e_i;

end