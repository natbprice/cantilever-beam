function [kopt,Ef,Epf,pre,Ef_var,Epf_var] = Create_Tradeoff_Curve(opts)
% Create tradeoff curve for expected performance vs probability of redesign
%
% Inputs:
% opts = structure with inputs and options
%
% Outputs:
% nopt = matrix of optimum designs
% Ef = column vector of expected cost function 
% Epf = column vector of expected probability of failure
% pre = column vector of probability of redesign
%==========================================================================

%% Unpack variables
method=opts.margins.method;
lb=opts.margins.lb;
ub=opts.margins.ub;
Srgts=opts.margins.Srgts;

%% Unpack surrogate models
model_f=Srgts.Ef{1};
xi_f=Srgts.Ef{2};
zi_f=Srgts.Ef{3};
kreq_f=Srgts.Ef{4};

% model_b=Srgts.beta{1};
% xi_b=Srgts.beta{2};
% zi_b=Srgts.beta{3};
% kreq_b=Srgts.beta{4};

%% Create vector of limits on probability of redesign
plim=linspace(0,0.5,6);
plim(1)=0.003;
n=length(plim);

%% Preallocate arrays
kopt=zeros(n,4);
Ef=zeros(n,1);
Epf=zeros(n,1);
Ef_var=zeros(n,1);
Epf_var=zeros(n,1);
pre=zeros(n,1);

%% Loop through limit values to create tradeoff curve (Pareto Front)
for i=1:n
    opts.margins.pre_bar=plim(i);
    
    % Call optimizer to find optimum safety factors
    [xmin,fmin,counteval,stopflag,out,bestever] = Safety_Margin_Opt(opts);
    xbest=bestever.x';
       
    % Fix klb and kub at bounds if not used
    switch lower(method)
        case {'mixed'}
        case {'safety'}
            xbest=[xbest(1),xbest(2),0.95,xbest(3)];
        case {'performance'}
            xbest=[xbest(1),0.95,xbest(2),xbest(3)];
    end
    kopt(i,:)=lb+(ub-lb).*xbest;
       
    % Surrogate prediction of objective function
    z=stk_predict(model_f,{xi_f,kreq_f},zi_f,xbest);
    Ef(i,1)=z.mean;
    Ef_var(i,1)=z.var;

%     % Surrogate prediction of reliability constraint
%     z=stk_predict(model_b,{xi_b,kreq_b},zi_b,xbest);
%     Epf(i,1)=normcdf(-z.mean);
%     Epf_var(i,1)=z.var;

    % Calculate probability of redesign
    pre(i,1)=normcdf(-kopt(i,2))+(1-normcdf(kopt(i,3)));   
    
    % Calculate probability of negative safety margin
    pd=makedist('Normal','mu',0,'sigma',1);
    pd_ini=truncate(pd,-kopt(i,2),kopt(i,3));
    Pf_n_ini=cdf(pd_ini,-kopt(i,1));
    Pf_n_re=normcdf(-kopt(i,4));
    Epf(i,1)=(1-pre(i,1))*Pf_n_ini+pre(i,1)*Pf_n_re;
    Epf_var(i,1)=0;
   
    % Display results
    fprintf('\n\nPre   Ewe   Ef\n')
    fprintf('%4.2f  %7.2e  %7.2e\n\n',[pre(i,1),Ef(i,1),Epf(i,1)])

end

%% Plot tradeoff curve
figure(1)
hold on
errorbar(pre,Ef,1.96*sqrt(Ef_var),'-o')
xlabel('Probability of redesign')
ylabel('Expected objective function value')

end

