function [dat] = Design_Process_Sim_Adapt(k,E,opts)
% Perform MCS of deterministic design process with adaptive sample size
%
% Inputs:
% k = standard deviation offsets vector
% E = structure with error model
% opts = structure with inputs and options
%
% Outputs:
% dat = structure with simulation data
%==========================================================================
%% Start function clock
tic;

%% Unpack standard deviation offset vector
kini=k(1);
klb=k(2);
kub=k(3);
kre=k(4);

%% Unpack variables
udet=opts.design.udet;
dim_x=opts.design.dim;
dim_u=opts.reliability.dim;
mstart=opts.error.mstart;
mmax=opts.error.mmax;
cov_target=opts.error.cov_target;
createHists=opts.fig.createHists;
calcPf=opts.reliability.calcPf;
DoE=opts.error.DoE_for_cond_sims;
pf_star=opts.margins.pf_bar;
E.returnKstd=false;

%% Initial design optimization
E.returnCalib=false;
E.i=0;
[xini,~]=Design_Opt(kini,udet,E,opts);

%% Map standard deviation offset (kini) to safety margin (nini)
nini=g_H(xini,udet,E,opts);

%% Calculate probability of redesign
pre=normcdf(-klb)+(1-normcdf(kub));

%% Calculate probability of negative safety margin
pd=makedist('Normal','mu',0,'sigma',1);
pd_ini=truncate(pd,-klb,kub);
Pf_n_ini=cdf(pd_ini,-kini);
Pf_n_re=normcdf(-kre);
Pf_n=(1-pre)*Pf_n_ini+pre*Pf_n_re;

%% Map bounds on standard deviation (klb,kub) to bounds on safety margin (nlb,nub)
i=0;                                % Use i=0 for initial DoE
model=E.model;                      % STK error model
Zi=E.Zi{i+1};                       % DoE
Ei=E.Ei{i+1};                       % Evaluations at DoE
kreq=E.kreq{i+1};                   % STK structure for fast evalutations
z=[xini,udet];

% Evaluate error model for initial design
epred=stk_predict(model,{Zi,kreq},Ei,z);
nlb_E=epred.mean-klb*sqrt(epred.var);
nub_E=epred.mean+kub*sqrt(epred.var);

% Evaluate constraint function of error bounds
if isinf(nlb_E)
    nlb=inf;
else
    nlb=g_H_of_E(xini,udet,nlb_E,opts);
end
if isinf(nub_E)
    nub=inf;
else
    nub=g_H_of_E(xini,udet,nub_E,opts);
end

%% Loop through GP trajectories
% Print header to screen and initialize figure
if opts.error.display
    fprintf('\n%7s %9s %12s %12s %7s %7s %7s %8s\n','i','t (min)','Ef_cov (%)','Pf_cov (%)','m1','m2','m','t pred')
    h=animatedline('Marker','o','LineStyle','-');
end

% Preallocate arrays
Pfini=zeros(mmax,1);
Bini=zeros(mmax,1);
Umpp_ini=zeros(mmax,dim_u);
GH_ini=zeros(mmax,1);
Nmeas_1=zeros(mmax,1);
Q=zeros(mmax,1);
Xre=zeros(mmax,dim_x);
Fini=zeros(mmax,1);
Fre=zeros(mmax,1);
Fini2=zeros(mmax,1);
Fre2=zeros(mmax,1);
GH_re=zeros(mmax,1);
Nmeas_2=zeros(mmax,1);
Xfinal=zeros(mmax,dim_x);
Pffinal=zeros(mmax,1);
Bfinal=zeros(mmax,1);
Umpp_final=zeros(mmax,dim_u);
Ffinal=zeros(mmax,1);
Ffinal2=zeros(mmax,1);
Nfinal=zeros(mmax,1);
nre=zeros(mmax,1);

% Initialize variables
Epf_cov=inf;
Ef_cov=inf;
if calcPf
    Pf_alpha_cov=inf;
else
    Pf_alpha_cov=0;
end
Ein=E;
Ein.i=1;
m=1;
i_start=1;
i_finish=mstart;

while (Pf_alpha_cov>cov_target(1) || Ef_cov>cov_target(2))  && m<mmax
    
    %% Evaluate epistemic realizations in parallel
    parfor i=i_start:i_finish
        %% Sample epistemic model uncertainty
        E=cond_sim_DoE([xini,udet;DoE],Ein,opts);

        %% Simulate future test (initial)
        E.returnCalib=false;
        GH_ini(i,1)=g_H(xini,udet,E,opts);

        %% Calculate safety margin (initial)
        Nmeas_1(i,1)=GH_ini(i,1)-0;

        %% Make redesign decision
        Q(i,1)=max(Nmeas_1(i,1)<nlb,Nmeas_1(i,1)>nub);

        if Q(i,1)==1
            %% Redesign design optimization
            E.returnCalib=true;
            [Xre(i,:),~]=Design_Opt(kre,udet,E,opts,xini);

            %% Calculate mean objective function value after redesign
            E.returnCalib=false;
            Fre(i,1)=f(Xre(i,:),udet,E,opts);
%             [Fre(i,1),Fre2(i,1)]=Objective_Function_MCS(Xre(i,:),E,opts);

            %% Simulate future test (redesign)
            E.returnCalib=false;
            GH_re(i,1)=g_H(Xre(i,:),udet,E,opts);

            %% Calculate safety margin (redesign)
            Nmeas_2(i,1)=GH_re(i,1)-0;
            
            %% Map standard deviation offset (kre) to safety margin (nre)
            E.returnCalib=true;
            nre(i,1)=g_H(Xre(i,:),udet,E,opts);
        end

        %% Calculate initial mean objective function value
        E.returnCalib=false;
        Fini(i,1)=f(xini,udet,E,opts);
%         [Fini(i,1),Fini2(i,1)]=Objective_Function_MCS(xini,E,opts);

        %% Calculate initial probability of failure
        if createHists
            if calcPf
                E.returnCalib=false;
                [Pfini(i,1),Bini(i,1),Umpp_ini(i,:)]=Reliability_Analysis(xini,E,opts);
            end
        end

        %% Calculate final design
        Xfinal(i,:)=(1-Q(i,1)).*xini+Q(i,1).*Xre(i,:);

        %% Calculate final reliability
        if calcPf
            E.returnCalib=false;
            [Pffinal(i,1),Bfinal(i,1),Umpp_final(i,:)]=Reliability_Analysis(Xfinal(i,:),E,opts);
        end

        %% Calculate final objective function value
        Ffinal(i,1)=(1-Q(i,1)).*Fini(i,1)+Q(i,1).*Fre(i,1);
        Ffinal2(i,1)=(1-Q(i,1)).*Fini2(i,1)+Q(i,1).*Fre2(i,1);

        %% Calculate final safety margin
        Nfinal(i,1)=(1-Q(i,1)).*Nmeas_1(i,1)+Q(i,1).*Nmeas_2(i,1);
    
    end
    
    %% Current number of epistemic samples
    m=i_finish;
    
    %% Define temporary arrays for calculations
    Q_temp=Q(1:m,1);
    Ffinal_temp=Ffinal(1:m,1);
    Ffinal_temp2=Ffinal2(1:m,1);
    Pffinal_temp=Pffinal(1:m,1);
    Xfinal_temp=norm_design(Xfinal(1:m,:),'x',0,opts);
    
    %% Calculate probability of violating probability of failure constraint
    Pf_alpha=sum(Pffinal>pf_star)/m;
    std_Pf_alpha=sqrt((Pf_alpha*(1-Pf_alpha))/m);
    
    %% Calculate mean probability of failure and mean objective function
    if sum(Q_temp)>=1 && sum(Q_temp)<m
        Ef=(1-pre)*mean(Ffinal_temp(Q_temp==0))+pre*mean(Ffinal_temp(Q_temp==1));
        Ef2=(1-pre)*mean(Ffinal_temp2(Q_temp==0))+pre*mean(Ffinal_temp2(Q_temp==1));
        Epf=(1-pre)*mean(Pffinal_temp(Q_temp==0))+pre*mean(Pffinal_temp(Q_temp==1));
        Ex=(1-pre)*mean(Xfinal_temp(Q_temp==0,:))+pre*mean(Xfinal_temp(Q_temp==1,:));
        
        m1=sum(Q_temp==0);
        m2=sum(Q_temp==1);
        
        var_f1=var(Ffinal_temp(Q_temp==0));
        var_f2=var(Ffinal_temp(Q_temp==1));
        std_Ef=sqrt((1-pre)^2*(var_f1/m1)+pre^2*var_f2/m2);
        
        var_Pf1=var(Pffinal_temp(Q_temp==0));
        var_Pf2=var(Pffinal_temp(Q_temp==1));
        std_Epf=sqrt((1-pre)^2*(var_Pf1/m1)+pre^2*var_Pf2/m2);
        
        var_f21=var(Ffinal_temp2(Q_temp==0));
        var_f22=var(Ffinal_temp2(Q_temp==1));
        std_Ef2=sqrt((1-pre)^2*(var_f21/m1)+pre^2*var_f22/m2);
        
        var_x1=var(Xfinal_temp(Q_temp==0,1));
        var_x2=var(Xfinal_temp(Q_temp==1,1));
        std_Ex1=sqrt((1-pre)^2*(var_x1/m1)+pre^2*var_x2/m2);
    elseif sum(Q_temp)==m
        % 100 percent redesign
        Ef=pre*mean(Ffinal_temp(Q_temp==1));
        Ef2=pre*mean(Ffinal_temp2(Q_temp==1));
        Epf=pre*mean(Pffinal_temp(Q_temp==1));
        Ex=pre*mean(Xfinal_temp(Q_temp==1,:));
        
        m1=0;
        m2=sum(Q_temp==1);
        
        var_f1=0;
        var_f2=var(Ffinal_temp(Q_temp==1));
        std_Ef=sqrt(pre^2*var_f2/m2);
        
        var_Pf1=0;
        var_Pf2=var(Pffinal_temp(Q_temp==1));
        std_Epf=sqrt(pre^2*var_Pf2/m2);
        
        var_f21=0;
        var_f22=var(Ffinal_temp2(Q_temp==1));
        std_Ef2=sqrt(pre^2*var_f22/m2);
        
        var_x1=0;
        var_x2=var(Xfinal_temp(Q_temp==1,1));
        std_Ex1=sqrt(pre^2*var_x2/m2);
    else
        % 0 percent redesign
        Ef=(1-pre)*mean(Ffinal_temp(Q_temp==0));
        Ef2=(1-pre)*mean(Ffinal_temp2(Q_temp==0));
        Epf=(1-pre)*mean(Pffinal_temp(Q_temp==0));
        Ex=(1-pre)*mean(Xfinal_temp(Q_temp==0,:));
                
        m1=sum(Q_temp==0);
        m2=0;
        
        var_f1=var(Ffinal_temp(Q_temp==0));
        var_f2=0;
        std_Ef=sqrt((1-pre)^2*(var_f1/m1));
        
        var_Pf1=var(Pffinal_temp(Q_temp==0));
        var_Pf2=0;
        std_Epf=sqrt((1-pre)^2*(var_Pf1/m1));
        
        var_f21=var(Ffinal_temp2(Q_temp==0));
        var_f22=0;
        std_Ef2=sqrt((1-pre)^2*(var_f21/m1));
        
        var_x1=var(Xfinal_temp(Q_temp==0,1));
        var_x2=0;
        std_Ex1=sqrt((1-pre)^2*(var_x1/m1));
    end

    %% Calculate CoV of expectations
    if m>=mstart
        Ef_cov=(std_Ef/Ef)*100;
        Ef2_cov=(std_Ef2/Ef2)*100;
        Ex1_cov=(std_Ex1/Ex(:,1))*100;
        Epf_cov=(std_Epf/Epf)*100;
        Pf_alpha_cov=(std_Pf_alpha/Pf_alpha)*100;
    end
    
    %% Restart epistemic sample loop
    i_start=i_finish+1;
    m1_pred=floor((sqrt((1-pre)*var_f1+pre*var_f2)/(Ef*0.01*cov_target(2)))^2);
%     m2_pred=floor((sqrt((1-pre)*var_Pf1+pre*var_Pf2)/(Epf*0.01*cov_target(1)))^2);
    m2_pred=floor((Pf_alpha*(1-Pf_alpha))/(Pf_alpha*0.01*cov_target(1))^2);
    m_pred=max(m1_pred,m2_pred);
    m_min=i_start+31;
    i_finish=min(max(m_pred,m_min),mmax);
    
    %% Print current iteration to screen
    if opts.error.display
        t0=toc/60;
        fprintf('%7.0i %9.1f %12.2f %12.2f %7.0i %7.0i %7.0i %8.1f\n',...
            [m,t0,(std_Ef/Ef)*100,(std_Pf_alpha/Pf_alpha)*100,m1_pred,m2_pred,i_finish,(t0/m)*i_finish])
        addpoints(h,m,(std_Pf_alpha/Pf_alpha)*100)
        drawnow
    end
end

%% Confidence intervals
Ef_95CI=[Ef-1.96*std_Ef,Ef+1.96*std_Ef];
Pf_alpha_95CI=[Pf_alpha-1.96*std_Pf_alpha,Pf_alpha+1.96*std_Pf_alpha];

%% Resample arrays
Pfini=Pfini(1:m,1);
Bini=Bini(1:m,1);
Umpp_ini=Umpp_ini(1:m,:);
Fini=Fini(1:m,1);
Fini2=Fini2(1:m,1);
GH_ini=GH_ini(1:m,1);
Nmeas_1=Nmeas_1(1:m,1);
Q=Q(1:m,1);
Xre=Xre(1:m,:);
Fre=Fre(1:m,1);
Fre2=Fre2(1:m,1);
GH_re=GH_re(1:m,1);
Nmeas_2=Nmeas_2(1:m,1);
Xfinal=Xfinal(1:m,:);
Pffinal=Pffinal(1:m,1);
Bfinal=Bfinal(1:m,1);
Umpp_final=Umpp_final(1:m,:);
Ffinal=Ffinal(1:m,1);
Ffinal2=Ffinal2(1:m,1);
Nfinal=Nfinal(1:m,1);

%% Final time
t=toc;

%% Create data structure
dat=v2struct(nini,nlb,nub,nre,m,Pfini,Bini,Umpp_ini,Fini,Fini2,GH_ini,Nmeas_1,Q,Xre,Fre,Fre2,...
    GH_re,Nmeas_2,xini,Xfinal,Pffinal,Bfinal,Umpp_final,Ffinal,Ffinal2,Nfinal,...
    Ef,Epf,pre,Ef_cov,std_Ef,Ef_95CI,Epf_cov,m,Ef2,Ex,Ef2_cov,Ex1_cov,Pf_alpha,Pf_alpha_cov,...
    Pf_alpha_95CI,std_Pf_alpha,Pf_n,Pf_n_ini,Pf_n_re,t);

%% Histograms
if createHists
    PlotHists(dat,opts);
end

end