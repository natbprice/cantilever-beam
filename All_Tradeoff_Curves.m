function [Mixed,Safety,Perf] = All_Tradeoff_Curves(opts)
% Create multiple tradeoff curves: safety, performance, mixed
%
% Inputs:
% opts = structure with inputs and options
%
% Outputs:
% Mixed = structure of results from possibly mixed redesign
% Safety = structure of results from redesign for safety curve
% Perf = structure of results from redesign for performance curve
%==========================================================================

opts.margins.method='mixed';
[Mixed.nopt, Mixed.Ef, Mixed.Epf, Mixed.pre, Mixed.Ef_var, Mixed.Epf_var]=Create_Tradeoff_Curve(opts);
save('Mixed','Mixed')

opts.margins.method='safety';
[Safety.nopt, Safety.Ef, Safety.Epf, Safety.pre, Safety.Ef_var, Safety.Epf_var]=Create_Tradeoff_Curve(opts);
save('Safety','Safety')

opts.margins.method='performance';
[Perf.nopt, Perf.Ef, Perf.Epf, Perf.pre, Perf.Ef_var, Perf.Epf_var]=Create_Tradeoff_Curve(opts);
save('Perf','Perf')

figure()
hold on
errorbar(Mixed.pre,Mixed.Ef,1.96*sqrt(Mixed.Ef_var),'-b')
errorbar(Safety.pre,Safety.Ef,1.96*sqrt(Safety.Ef_var),'-r')
errorbar(Perf.pre,Perf.Ef,1.96*sqrt(Perf.Ef_var),'-g')
xlabel('Probability of redesign')
ylabel('Expected cross sectional area (in^2)')
legend('Mixed','Safety','Performance')

end

