function [Mixed,Safety,Perf] = MCS_tradeoff_curves(E,opts)
% Perform MCS for each point on tradeoff curve
%
% Inputs:
% E = structure with error model
% opts = structure with inputs and options
%
% Outputs:
% Mixed = updated structure with MCS results
% Safety = updated structure with MCS results
% Perf = updated structure with MCS results
%==========================================================================

%% Load tradeoff curve based on surrogates
S=load([pwd,'\results\','tradeoff_curves']);
Mixed=S.Mixed;
Perf=S.Perf;
Safety=S.Safety;

%% Header
fprintf('%6s %6s %9s %10s %8s %10s %14s %6s\n','i','m','t (min)','Ef','Ef_cov','Pf_alpha','Pf_alpha_cov','pre');

%% Mixed redesign
n=Mixed.nopt;
for i=1:6
       
    % Evaluate model at ith point in DoE
    [dat{i}] = Design_Process_Sim_Adapt(n(i,:),E,opts);
    v2struct(dat{i});
    
    % Display results for ith point in DoE
    fprintf('%6.0i %6.0i %9.1f %10.2e %8.2f %10.3f %14.2f %6.2f\n',[i,m,t/60,Ef,Ef_cov,Pf_alpha,Pf_alpha_cov,pre])
end
Mixed.MCS=dat;

%% Redesign for safety
n=Safety.nopt;
for i=1:6
       
    % Evaluate model at ith point in DoE
    [dat{i}] = Design_Process_Sim_Adapt(n(i,:),E,opts);
    v2struct(dat{i});
    
    % Display results for ith point in DoE
    fprintf('%6.0i %6.0i %9.1f %10.2e %8.2f %10.3f %14.2f %6.2f\n',[i,m,t/60,Ef,Ef_cov,Pf_alpha,Pf_alpha_cov,pre])
end
Safety.MCS=dat;

%% Redesign for performance
n=Perf.nopt;
for i=1:6
       
    % Evaluate model at ith point in DoE
    [dat{i}] = Design_Process_Sim_Adapt(n(i,:),E,opts);
    v2struct(dat{i});
    
    % Display results for ith point in DoE
    fprintf('%6.0i %6.0i %9.1f %10.2e %8.2f %10.3f %14.2f %6.2f\n',[i,m,t/60,Ef,Ef_cov,Pf_alpha,Pf_alpha_cov,pre])
end
Perf.MCS=dat;

%% Plots
close all

figure(1)
hold on
data=Mixed;
for i=1:6
    pre(i)=data.MCS{i}.pre;
    Pf_alpha(i)=data.MCS{i}.Pf_alpha;
    std_Pf_alpha(i)=data.MCS{i}.std_Pf_alpha;
end
errorbar(pre,Pf_alpha,1.96*std_Pf_alpha,'-bo','MarkerFaceColor','b','MarkerSize',1);
xlim([-0.01,0.51])
xlabel('Probability of redesign')
ylabel('Probability of constraint violation')
ylim([0,0.15])

figure(2)
hold on
data=Safety;
for i=1:6
    pre(i)=data.MCS{i}.pre;
    Pf_alpha(i)=data.MCS{i}.Pf_alpha;
    std_Pf_alpha(i)=data.MCS{i}.std_Pf_alpha;
end
errorbar(pre,Pf_alpha,1.96*std_Pf_alpha,'-bo','MarkerFaceColor','b','MarkerSize',1)
xlim([-0.01,0.51])
xlabel('Probability of redesign')
ylabel('Probability of constraint violation')
ylim([0,0.15])

figure(3)
hold on
data=Perf;
for i=1:6
    pre(i)=data.MCS{i}.pre;
    Pf_alpha(i)=data.MCS{i}.Pf_alpha;
    std_Pf_alpha(i)=data.MCS{i}.std_Pf_alpha;
end
errorbar(pre,Pf_alpha,1.96*std_Pf_alpha,'-bo','MarkerFaceColor','b','MarkerSize',1)
xlim([-0.01,0.51])
xlabel('Probability of redesign')
ylabel('Probability of constraint violation')
ylim([0,0.15])


figure(1)
hold on
data=Mixed;
for i=1:6
    pre(i)=data.MCS{i}.pre;
    Pf_n_mcs(i)=sum(data.MCS{i}.Nfinal<0)/data.MCS{i}.m;
    std_Pf_n_mcs(i)=sqrt((Pf_n_mcs(i)*(1-Pf_n_mcs(i)))/data.MCS{i}.m);
end
errorbar(pre,Pf_n_mcs,1.96*std_Pf_n_mcs,'-go','MarkerFaceColor','g','MarkerSize',1);
plot(xlim,[0.05,0.05],'--k')
legend('P_E[P_f(X_{final})> p_f*]','P_E[G(X_{final},u_{det})< 0]','\alpha')

figure(2)
hold on
data=Safety;
for i=1:6
    pre(i)=data.MCS{i}.pre;
    Pf_n_mcs(i)=sum(data.MCS{i}.Nfinal<0)/data.MCS{i}.m;
    std_Pf_n_mcs(i)=sqrt((Pf_n_mcs(i)*(1-Pf_n_mcs(i)))/data.MCS{i}.m);
end
errorbar(pre,Pf_n_mcs,1.96*std_Pf_n_mcs,'-go','MarkerFaceColor','g','MarkerSize',1);
plot(xlim,[0.05,0.05],'--k')
legend('P_E[P_f(X_{final})> p_f*]','P_E[G(X_{final},u_{det})< 0]','\alpha')

figure(3)
hold on
data=Perf;
for i=1:6
    pre(i)=data.MCS{i}.pre;
    Pf_n_mcs(i)=sum(data.MCS{i}.Nfinal<0)/data.MCS{i}.m;
    std_Pf_n_mcs(i)=sqrt((Pf_n_mcs(i)*(1-Pf_n_mcs(i)))/data.MCS{i}.m);
end
errorbar(pre,Pf_n_mcs,1.96*std_Pf_n_mcs,'-go','MarkerFaceColor','g','MarkerSize',1);
plot(xlim,[0.05,0.05],'--k')
legend('P_E[P_f(X_{final})> p_f*]','P_E[G(X_{final},u_{det})< 0]','\alpha')

end