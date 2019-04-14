function [g_U,g_E,Var_U,Var_E,R] = VarRatio(E,opts,x)
% Calculate ratio of epistemic to aleatory variance
%
% Inputs:
% x = vector of design variables (optional)
%
% Outputs:
% g_U = G(x,U) | G(.,.)=\bar{g}(.,.)
% g_E = G(x,U) | U=\bar{u}
% Var_U = Var(g_U)
% Var_E = Var(g_E)
% R = Var_E / Var_U
%==========================================================================

%% Number of samples
n=1e2;      % Aleatory
m=1e2;      % Epistemic

%% Random design point
if nargin==2
    dim_x=opts.design.dim;
    x=rand(1,dim_x);
end

%% Calculate aleatory variance conditional on mean epistemic uncertainty
E.returnCalib=false;
dim_u=opts.reliability.dim;
U_norm=myTinv(normrnd(0,1,n,dim_u),opts);
g_U=zeros(n,1);
for j=1:n
    g_U(j,1)=g_H(x,U_norm(j,:),E,opts);
end
Var_U=var(g_U);

%% Calculate epistemic variance conditional on mean aleatory uncertainty
Ein=E;
Ein.i=1;
g_E=zeros(m,1);
DoE=opts.error.DoE_for_cond_sims;
parfor i=1:m
    E=cond_sim_DoE(DoE,Ein,opts);
    E.returnCalib=false;
    g_E(i,1)=g_H(x,myTinv([0,0],opts),E,opts);
end
Var_E=var(g_E);

%% Variance ratio
R=Var_E/Var_U;

% %% Plot overlapping histograms
% figure(1)
% clf
% hold on
% fc1='b';
% fc2='y';
% ec='k';
% falph=0.4;
% nbins=20;
% x_edge=linspace(min([g_U;g_E]),max([g_U;g_E]),nbins+1);
% x_edge=[x_edge(1:nbins),inf];
% y=histc(g_U,x_edge);
% h1=bar(x_edge(1:nbins),y(1:nbins),'histc');
% set(h1,'facecolor',fc1,'EdgeColor',ec,'linestyle','-')
% y=histc(g_E,x_edge);
% h2=bar(x_edge(1:nbins),y(1:nbins),'histc');
% set(h1,'facecolor',fc2,'EdgeColor',ec,'linestyle','-')
% set(findobj(gca,'type','patch'),'facea',falph)

end

