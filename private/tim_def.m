function d=tim_def(x_norm,u_norm,opts)
% Timoshenko model of tip deflection for cantilever beam
%
% Inputs:
% x = vector of design variables
% u = vector of aleatory variables
% opts = structure with inputs and options
%
% Outputs:
% d = tip deflection
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

% Eq. xiv, pg. 60, Aircraft Structures for Engineering Students, 4th
% ed, Megson, 2007
dx=((Px*x.^3)./(6*E*Iz)) - ((Px*L^2*x)./(2*E*Iz)) + ((Px*L^3)./(3*E*Iz)) + ((Px.*h.^2)./(8*Iz*G))*(L-x);
dy=((Py*x.^3)./(6*E*Iy)) - ((Py*L^2*x)./(2*E*Iy)) + ((Py*L^3)./(3*E*Iy)) + ((Py.*b.^2)./(8*Iy*G))*(L-x);
d=sqrt(dx.^2+dy.^2);
end

