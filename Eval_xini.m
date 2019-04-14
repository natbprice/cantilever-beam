function [data] = Eval_xini(data,opts,E)
% Evaluate true probability of failure and safety margin for tradeoff curve

for j=1:6
    n=data.nopt(j,:);
    [xini(j,:),~]=Design_Opt(n(1),opts.design.udet,E,opts);
    [Pfini(j,1),Bini(j,1)]=Reliability_Analysis(xini(j,:),E,opts,'high');
    [~,~,~,nHini(j,1)]=Reliability_Analysis(xini(j,:),E,opts,'high');
    
    %% Calculate probability of redesign
    i=0;                                % Use i=0 for initial DoE
    model=E.model;                      % STK error model
    Zi=E.Zi{i+1};                       % DoE
    Ei=E.Ei{i+1};                       % Evaluations at DoE
    kreq=E.kreq{i+1};                   % STK structure for fast evalutations
    z=[xini(j,:),opts.design.udet];

    epred=stk_predict(model,{Zi,kreq},Ei,z);
    nlb_E=epred.mean-n(2)*sqrt(epred.var);
    nub_E=epred.mean+n(3)*sqrt(epred.var);
    
    if isinf(nlb_E)
        nlb(j,1)=inf;
    else
        nlb(j,1)=g_H_of_E(xini(j,:),opts.design.udet,nlb_E,opts);
    end
    if isinf(nub_E)
        nub(j,1)=inf;
    else
        nub(j,1)=g_H_of_E(xini(j,:),opts.design.udet,nub_E,opts);
    end

    Q(j,1)=max(nHini(j,1)<nlb(j,1),nHini(j,1)>nub(j,1));
    
end

data.xini=xini;
data.Pf=Pfini;
data.B=Bini;
data.n=[nHini,nlb,nub];
data.Q=Q;
