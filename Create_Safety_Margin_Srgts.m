function [srgts,e_pf,e_f]=Create_Safety_Margin_Srgts(method)
% Create surrogates as a function of safety margins
%
% Inputs:
% method = redesign method
%
% Outputs:
% srgts = structure with surrogate models
% e_pf = error in mean probability of failure
% e_f = error in expected cost
%==========================================================================

close all

%% Switch case for different redesign methods
outDir=[pwd,'\safety_margin_doe\'];
switch lower(method)
    case {'mixed'}
%         dim=4;
%         S=load([outDir,'margins_DoE_',method],'Ef','Epf','n_norm');
%         Ef1=S.Ef;
%         ind=length(Ef1);
%         Epf1=S.Epf;
%         n_norm1=S.n_norm(1:ind,:);
%         S=load([outDir,'margins_DoE_','safety'],'Ef','Epf','n_norm');
%         Ef2=S.Ef;
%         ind=length(Ef2);
%         Epf2=S.Epf;
%         n_norm2=S.n_norm(1:ind,:);
%         n_norm2=[n_norm2(:,1),n_norm2(:,2),ones(ind,1),n_norm2(:,3)];
%         S=load([outDir,'margins_DoE_','performance'],'Ef','Epf','n_norm');
%         Ef3=S.Ef;
%         ind=length(Ef3);
%         Epf3=S.Epf;
%         n_norm3=S.n_norm(1:ind,:);
%         n_norm3=[n_norm3(:,1),ones(ind,1),n_norm3(:,2),n_norm3(:,3)];
%         Ef=[Ef1;Ef2;Ef3];
%         Epf=[Epf1;Epf2;Epf3];
%         n_norm=[n_norm1;n_norm2;n_norm3];

        dim=4;
        S=load([outDir,'margins_DoE_',method],'Ef','Pf_n','n_norm','Pf_n_cov');
        Ef=S.Ef;
        Pf_n_cov=S.Pf_n_cov/100;
        ind=length(Ef);
        Pf_n=S.Pf_n;
        n_norm=S.n_norm(1:ind,:);
        
%         warning('Removing points with Pf=0')
%         ind=Pf_n~=0;
%         Pf_n=Pf_n(ind,:);
%         n_norm=n_norm(ind,:);
%         Ef=Ef(ind,:);
%         Pf_n_cov=Pf_n_cov(ind,:);
    case {'safety'}
        dim=3;
        S=load([outDir,'margins_DoE_',method],'Ef','Epf','n_norm','Ef_cov','Epf_cov');
        Ef=S.Ef;
        Ef_cov=S.Ef_cov/100;
        Epf_cov=S.Epf_cov/100;
        ind=length(Ef);
        Epf=S.Epf;
        n_norm=S.n_norm(1:ind,:);
    case {'performance'}
        dim=3;
        S=load([outDir,'margins_DoE_',method],'Ef','Epf','n_norm','Ef_cov','Epf_cov');
        Ef=S.Ef;
        Ef_cov=S.Ef_cov/100;
        Epf_cov=S.Epf_cov/100;
        ind=length(Ef);
        Epf=S.Epf;
        n_norm=S.n_norm(1:ind,:);
    otherwise
        error('Method not recognized')
end

%% Fit reliability index model
% Set covariance function and model order
model=stk_model('stk_gausscov_aniso',dim);
model.lognoisevariance=NaN;

% Compute an initial guess for the covariance parameters and a reasonable
% log-variance for a small "regularization noise"
box = [zeros(1,dim);ones(1,dim)];
[param0, model.lognoisevariance]=stk_param_init(model,n_norm,-norminv(Pf_n),box);

% Estimate covariance parameters from data
model.param=stk_param_estim(model,n_norm,-norminv(Pf_n),param0);

% Output models
srgts.beta{1}=model;
srgts.beta{2}=n_norm;
srgts.beta{3}=-norminv(Pf_n);
srgts.beta{4}=stk_kreq_qr(model,n_norm);

% Perform leave-one-out cross validation
if nargout==3
    parfor i=17:length(n_norm)
        ind=true(length(n_norm),1);
        ind(i)=false;
    %     model.param=stk_param_estim(model,n_norm(ind,:),-norminv(Epf(ind,:)),param0);
        z=stk_predict(model,n_norm(ind,:),-norminv(Pf_n(ind,:)),n_norm(i,:));
        e_pf_rel(i,1)=(Pf_n(i,:)-normcdf(-z.mean))/Pf_n(i,:);
        e_pf_abs(i,1)=(Pf_n(i,:)-normcdf(-z.mean));
    end
    e_pf=[e_pf_rel,e_pf_abs];
    e_pf_90_CI=[quantile(e_pf_rel,0.05),quantile(e_pf_rel,0.95)]
    Pf_cov_90_CI=[quantile(Pf_n_cov,0.05),quantile(Pf_n_cov,0.95)]

    figure()
    hold on
    hist(e_pf_rel,20)
    plot(repmat(quantile(e_pf_rel,0.05),2,1),ylim,'--r')
    h=plot(repmat(quantile(e_pf_rel,0.95),2,1),ylim,'--r');
    xlabel('Relative error')
    title('Mean probability of failure')
    ylabel('Number of occurences')
    legend(h,'90% CI','location','southeast')
    
    figure()
    hold on
    ecdf(e_pf_abs)
    xlim([quantile(e_pf_abs,0.025),quantile(e_pf_abs,0.975)])
    plot(repmat(quantile(e_pf_abs,0.05),2,1),ylim,'--r')
    h=plot(repmat(quantile(e_pf_abs,0.95),2,1),ylim,'--r');
    xlabel('Absolute error')
    title('Mean probability of failure')
    ylabel('CDF')
    legend(h,'90% CI','location','southeast')
    
    figure()
    hold on
    ecdf(e_pf_rel)
    xlim([quantile(e_pf_rel,0.025),quantile(e_pf_rel,0.975)])
    plot(repmat(quantile(e_pf_rel,0.05),2,1),ylim,'--r')
    h=plot(repmat(quantile(e_pf_rel,0.95),2,1),ylim,'--r');
    xlabel('Relative error')
    title('Mean probability of failure')
    ylabel('CDF')
    legend(h,'90% CI','location','southeast')
    
%     figure()
%     hold on
%     ecdf(Pf_n_cov)
%     xlim([quantile(Pf_n_cov,0.025),quantile(Pf_n_cov,0.975)])
%     plot(repmat(quantile(Pf_n_cov,0.05),2,1),ylim,'--r')
%     h=plot(repmat(quantile(Pf_n_cov,0.95),2,1),ylim,'--r');
%     xlabel('Coefficient of variation')
%     title('Mean probability of failure')
%     ylabel('CDF')
%     legend(h,'90% CI','location','southeast')

%     figure()
%     hold on
%     hist(e_pf_abs,20)
%     plot(repmat(quantile(e_pf_abs,0.025),2,1),ylim,'--r')
%     h=plot(repmat(quantile(e_pf_abs,0.975),2,1),ylim,'--r');
%     xlabel('Mean probability of failure - Absolute error')
%     ylabel('Number of occurences')
%     legend(h,'95% CI')
end

% %% Fit probability of failure model
% % Set covariance function and model order
% model=stk_model('stk_gausscov_aniso',dim);
% model.lognoisevariance=NaN;
% 
% % Compute an initial guess for the covariance parameters and a reasonable
% % log-variance for a small "regularization noise"
% box = [zeros(1,dim);ones(1,dim)];
% [param0, model.lognoisevariance]=stk_param_init(model,n_norm,Pf_n,box);
% 
% % Estimate covariance parameters from data
% model.param=stk_param_estim(model,n_norm,Pf_n,param0);
% 
% % Output models
% srgts.beta{1}=model;
% srgts.beta{2}=n_norm;
% srgts.beta{3}=Pf_n;
% srgts.beta{4}=stk_kreq_qr(model,n_norm);
% 
% % Perform leave-one-out cross validation
% if nargout==3
%     parfor i=17:length(n_norm)
%         ind=true(length(n_norm),1);
%         ind(i)=false;
%     %     model.param=stk_param_estim(model,n_norm(ind,:),-norminv(Epf(ind,:)),param0);
%         z=stk_predict(model,n_norm(ind,:),Pf_n(ind,:),n_norm(i,:));
%         e_pf_rel(i,1)=(Pf_n(i,:)-z.mean)/Pf_n(i,:);
%         e_pf_abs(i,1)=(Pf_n(i,:)-z.mean);
%     end
%     e_pf=[e_pf_rel,e_pf_abs];
%     e_pf_90_CI=[quantile(e_pf_rel,0.05),quantile(e_pf_rel,0.95)]
%     Pf_cov_90_CI=[quantile(Pf_n_cov,0.05),quantile(Pf_n_cov,0.95)]
% 
%     figure()
%     hold on
%     hist(e_pf_rel,20)
%     plot(repmat(quantile(e_pf_rel,0.05),2,1),ylim,'--r')
%     h=plot(repmat(quantile(e_pf_rel,0.95),2,1),ylim,'--r');
%     xlabel('Relative error')
%     title('Mean probability of failure')
%     ylabel('Number of occurences')
%     legend(h,'90% CI','location','southeast')
%     
%     figure()
%     hold on
%     ecdf(e_pf_abs)
%     xlim([quantile(e_pf_abs,0.025),quantile(e_pf_abs,0.975)])
%     plot(repmat(quantile(e_pf_abs,0.05),2,1),ylim,'--r')
%     h=plot(repmat(quantile(e_pf_abs,0.95),2,1),ylim,'--r');
%     xlabel('Absolute error')
%     title('Mean probability of failure')
%     ylabel('CDF')
%     legend(h,'90% CI','location','southeast')
%     
%     figure()
%     hold on
%     ecdf(e_pf_rel)
%     xlim([quantile(e_pf_rel,0.025),quantile(e_pf_rel,0.975)])
%     plot(repmat(quantile(e_pf_rel,0.05),2,1),ylim,'--r')
%     h=plot(repmat(quantile(e_pf_rel,0.95),2,1),ylim,'--r');
%     xlabel('Relative error')
%     title('Mean probability of failure')
%     ylabel('CDF')
%     legend(h,'90% CI','location','southeast')
%     
%     figure()
%     hold on
%     ecdf(Pf_n_cov)
%     xlim([quantile(Pf_n_cov,0.025),quantile(Pf_n_cov,0.975)])
%     plot(repmat(quantile(Pf_n_cov,0.05),2,1),ylim,'--r')
%     h=plot(repmat(quantile(Pf_n_cov,0.95),2,1),ylim,'--r');
%     xlabel('Coefficient of variation')
%     title('Mean probability of failure')
%     ylabel('CDF')
%     legend(h,'90% CI','location','southeast')
% 
% %     figure()
% %     hold on
% %     hist(e_pf_abs,20)
% %     plot(repmat(quantile(e_pf_abs,0.025),2,1),ylim,'--r')
% %     h=plot(repmat(quantile(e_pf_abs,0.975),2,1),ylim,'--r');
% %     xlabel('Mean probability of failure - Absolute error')
% %     ylabel('Number of occurences')
% %     legend(h,'95% CI')
% end

%% Fit Ef model
% Set covariance function and model order
model=stk_model('stk_gausscov_aniso',dim);
model.lognoisevariance=NaN;

% Compute an initial guess for the covariance parameters and a reasonable
% log-variance for a small "regularization noise"
box = [zeros(1,dim);ones(1,dim)];
[param0, model.lognoisevariance]=stk_param_init(model,n_norm,Ef,box);

% Estimate covariance parameters from data
model.param=stk_param_estim(model,n_norm,Ef,param0);

% Output models
srgts.Ef{1}=model;
srgts.Ef{2}=n_norm;
srgts.Ef{3}=Ef;
srgts.Ef{4}=stk_kreq_qr(model,n_norm);

% Perform leave-one-out cross validation
if nargout==3
    parfor i=17:length(n_norm)
        ind=true(length(n_norm),1);
        ind(i)=false;
    %     model.param=stk_param_estim(model,n_norm(ind,:),Ef(ind,:),param0);
        z=stk_predict(model,n_norm(ind,:),Ef(ind,:),n_norm(i,:));
        e_f_rel(i,1)=(Ef(i,:)-z.mean)/Ef(i,:);
        e_f_abs(i,1)=(Ef(i,:)-z.mean);
    end
    e_f=[e_f_rel,e_f_abs];
    e_f_90_CI=[quantile(e_f_rel,0.05),quantile(e_f_rel,0.95)]
%     Ef_cov_90_CI=[quantile(Ef_cov,0.05),quantile(Ef_cov,0.95)]

    figure()
    hold on
    hist(e_f_rel,20)
    plot(repmat(quantile(e_f_rel,0.05),2,1),ylim,'--r')
    h=plot(repmat(quantile(e_f_rel,0.95),2,1),ylim,'--r');
    xlabel('Relative error')
    title('Mean GLOW')
    ylabel('Number of occurences')
    legend(h,'90% CI','location','southeast')
    
    figure()
    hold on
    ecdf(e_f_rel)
    xlim([quantile(e_f_rel,0.025),quantile(e_f_rel,0.975)])
    plot(repmat(quantile(e_f_rel,0.05),2,1),ylim,'--r')
    h=plot(repmat(quantile(e_f_rel,0.95),2,1),ylim,'--r');
    xlabel('Relative error')
    title('Mean GLOW')
    ylabel('CDF')
    legend(h,'90% CI','location','southeast')
    
%     figure()
%     hold on
%     ecdf(Ef_cov)
%     xlim([quantile(Ef_cov,0.025),quantile(Ef_cov,0.975)])
%     plot(repmat(quantile(Ef_cov,0.05),2,1),ylim,'--r')
%     h=plot(repmat(quantile(Ef_cov,0.95),2,1),ylim,'--r');
%     xlabel('Coefficient of variation')
%     title('Mean GLOW')
%     ylabel('CDF')
%     legend(h,'90% CI','location','southeast')
    
%     figure()
%     hold on
%     hist(e_f_abs,20)
%     plot(repmat(quantile(e_f_abs,0.025),2,1),ylim,'--r')
%     h=plot(repmat(quantile(e_f_abs,0.975),2,1),ylim,'--r');
%     xlabel('Mean mass - Absolute error')
%     ylabel('Number of occurences')
%     legend(h,'95% CI')
end

end