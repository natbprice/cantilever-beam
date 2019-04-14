function [] = Safety_Margin_DoE(method,n,restart,newDoE,E,opts)
% Create DoE for fitting surrogates as a function of safety margins
%
% Inputs:
% method = redesign method
% n = number of points in DoE
% restart = flag to restart from existing saved data file
% E = structure with error model
% opts = structure with inputs and options
%
% Outputs:
% -- All variables are saved to data file
%==========================================================================

if restart
    i=[];
    load([pwd,'\safety_margin_doe\','margins_DoE_',method])
    outDir=[pwd,'\safety_margin_doe\'];
    istart=i+1;
else
    istart=1;
    
    %% Unpack variables
    xi=E.Zi{1};
    zi=E.Ei{1};
    model=E.model;

    %% Path to save data file
    outDir=[pwd,'\safety_margin_doe\'];

    %% Switch case for different redesign methods
    switch lower(method)
        case {'mixed'}
            dim=4;
            ind=logical([1,1,1,1]);
        case {'safety'}
            dim=3;
            ind=logical([1,1,0,1]);
        case {'performance'}
            dim=3;
            ind=logical([1,0,1,1]);
        otherwise
            error('Method not recognized')
    end
    
    if newDoE
        %% LHS Design
        n_norm_lhs=lhsdesign(n,dim,'criterion','maximin','iterations',5000);

        %% Corner points
        lb=ones(1,dim);
        ub=2*ones(1,dim);
        xff=fullfact(ub);
        n_norm_ff=bsxfun(@rdivide,bsxfun(@minus,xff,lb),(ub-lb));

        %% Combined DoE
        n_norm=[n_norm_ff;n_norm_lhs];

        %% Scale to design space
        lb=opts.margins.lb(ind);
        ub=opts.margins.ub(ind);
        n=bsxfun(@plus,lb,bsxfun(@times,n_norm,(ub-lb)));
    else
        %% Load existing DoE
        S=load([outDir,'n']);
        n=S.n;
        S=load([outDir,'n_norm']);
        n_norm=S.n_norm;
    end
    
    %% Pad DoE with inf values
    n(:,~ind)=inf;
    n_norm=n_norm(:,ind);
end

%% Evaluate model at points in DoE
fprintf('%6s %6s %9s %10s %8s %16s %6s\n','i','m','t (min)','Ef','Ef_cov','-norminv(Pf_n)','pre');
p=length(n);
for i=istart:p
        
    % Evaluate model at ith point in DoE
    [dat] = Design_Process_Sim_Adapt(n(i,:),E,opts);
    Ef(i,1)=dat.Ef;
    Pf_n(i,1)=dat.Pf_n;
    pre(i,1)=dat.pre;
    Ef_cov(i,1)=dat.Ef_cov;
    m(i,1)=dat.m;
    t(i,1)=dat.t;
    
    % Display results for ith point in DoE
    fprintf('%6.0i %6.0i %9.1f %10.2e %8.2f %16.2f %6.2f\n',[i,m(i,1),t(i,1)/60,Ef(i,1),Ef_cov(i,:),-norminv(Pf_n(i,1)),pre(i,1)])
    
    % Save all variables to data file
    save([outDir,'margins_DoE_',method])
end

end