function [f2] = RBDO_wrapper(x_norm,P)
% Calculate penalized objective function considering probability of failure
%
% Inputs:
% x_norm = normalized vector of design variables
%
% Outputs:
% f = penalized objective function
%==========================================================================

opts=P{1};
Ein=P{2};
model_fidelity=P{3};

w1=opts.margins.w1;
pf_bar=opts.margins.pf_bar;

x_norm=x_norm';
f2=zeros(1,length(x_norm(:,1)));

parfor i=1:length(x_norm(:,1))

    %% RBDO
    Ef=f(x_norm(i,:),[0,0],Ein,opts);
    [Pf,~,~]=Reliability_Analysis(x_norm(i,:),Ein,opts,model_fidelity);
    
    g1=Pf./pf_bar-1;
    f2(i)=Ef+w1*max(0,g1);
    
end

end

function g2=g(x,E,opts,model_fidelity)

    switch lower(model_fidelity)
        case {'low'}
            e_i=0;
            myFun='g_H_of_E';
            varIn={e_i,opts};
        case {'low_upd'}
            myFun='g_H';
            varIn={E,opts};
        case {'high'}
            myFun='g_H_true';
            varIn={opts};
        otherwise
            error('Unrecognized model fidelity')
    end
    
    [~,g2]=myEval(myFun,x,varIn);
end

function [g1,g2]=myEval(fun,x,varargin)
    varargin=varargin{1};
    nVarargs = length(varargin);
    if strcmp(fun,'g_H')
        [g1,g2]=feval(fun,x,0,varargin{1},varargin{2});
    else
        if nVarargs==1
            [~,~,g1,g2]=feval(fun,x,0,varargin{1});
        elseif nVarargs==2
            [~,~,g1,g2]=feval(fun,x,0,varargin{1},varargin{2});
        end
    end
end