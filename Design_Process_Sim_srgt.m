function [f] = Design_Process_Sim_srgt(x_norm,opts)
% Call STK surrogate models in place of Design_Process_Sim_Adapt.m
%
% Inputs:
% x_norm = normalized vector of standard deviation offsets
% opts = structure with inputs and options
%
% Outputs:
% f = penalized objective function
%==========================================================================

%% Unpack variables
Srgts=opts.margins.Srgts;
method=opts.margins.method;
lb=opts.margins.lb;
ub=opts.margins.ub;
alpha_bar=opts.margins.alpha_bar;
pre_bar=opts.margins.pre_bar;
w1=opts.margins.w1;
w2=opts.margins.w2;

%% Unpack surrogate models
model_f=Srgts.Ef{1};
xi_f=Srgts.Ef{2};
zi_f=Srgts.Ef{3};
kreq_f=Srgts.Ef{4};

% model_b=Srgts.beta{1};
% xi_b=Srgts.beta{2};
% zi_b=Srgts.beta{3};
% kreq_b=Srgts.beta{4};

%% Fix klb and kub at bounds if not used
switch lower(method)
    case {'mixed'}
    case {'safety'}
        x_norm=[x_norm(1,:);x_norm(2,:);0.95*ones(size(x_norm(1,:)));x_norm(3,:)];       
    case {'performance'}
        x_norm=[x_norm(1,:);0.95*ones(size(x_norm(1,:)));x_norm(2,:);x_norm(3,:)];
end

%% Surrogate prediction of objective function
z=stk_predict(model_f,{xi_f,kreq_f},zi_f,x_norm');
Ef=z.mean;

% %% Surrogate prediction of reliability constraint
% z=stk_predict(model_b,{xi_b,kreq_b},zi_b,x_norm');
% Pf=normcdf(-z.mean);

%% Calculate probability of redesign
n=bsxfun(@plus,lb,bsxfun(@times,x_norm',(ub-lb)));
pre=normcdf(-n(:,2))+(1-normcdf(n(:,3)));

%% Calculate probability of negative safety margin
pd=makedist('Normal','mu',0,'sigma',1);
[r,c]=size(n);
Pf=zeros(r,1);
parfor i=1:r
    pd_ini=truncate(pd,-n(i,2),n(i,3));
    Pf_n_ini=cdf(pd_ini,-n(i,1));
    Pf_n_re=normcdf(-n(i,4));
    Pf(i,1)=(1-pre(i,1))*Pf_n_ini+pre(i,1)*Pf_n_re;
end

%% Define constraint functions
g1=Pf./alpha_bar-1;
g2=pre./pre_bar-1;

%% Return penalized objective function
f=Ef+w1*max(0,g1)+w2*max(0,g2);
% f=Ef+w1*max(0,g1)+w2*abs(g2);
f=f';

end


