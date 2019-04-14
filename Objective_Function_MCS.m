function [EU_Mtot,EU_Minert] = Objective_Function_MCS(x,E,opts,model_fidelity)
% Calculate mean value of objective function w.r.t. aleatory uncertainty
%
% Inputs:
% x  = matrix of design variables
% E = structure with error model
% opts = structure with inputs and options
%
% Outputs:
% EU_f = mean objective function value
%==========================================================================

%% Optional model fidelity specification
if nargin==4
    switch lower(model_fidelity)
        case {'low'}
            e_i=0;
            myFun='g_H_of_E';
            varIn={x,e_i,opts};
        case {'low_upd'}
            myFun='f';
            varIn={x,E,opts};
        case {'high'}
            myFun='g_H_true';
            varIn={x,opts};
        otherwise
            error('Unrecognized model fidelity')
    end
else
    myFun='f';
    varIn={x,E,opts};
end

%% Unpack options
dim_u=opts.reliability.dim;
nmax=opts.obj.nmax;
nstart=opts.obj.nstart;
cov_target=opts.obj.cov;

%% Loop through aleatory realizations
U=zeros(nmax,dim_u);
M_tot=zeros(nmax,1);
M_inert=zeros(nmax,1);
i=1;
EU_Mtot_cov=inf;
while EU_Mtot_cov>cov_target && i<=nmax
    U(i,:)=normrnd(0,1,1,dim_u);
    [M_tot(i,1),M_inert(i,1)]=myEval(myFun,myTinv(U(i,:),opts),varIn);
%     [M_tot(i,1),M_inert(i,1)]=f(x,myTinv(U(i,:),opts),E,opts);
    
    %% Calculate outputs
    EU_Mtot=mean(M_tot(1:i,1));
    EU_Minert=mean(M_inert(1:i,1));

    %% Calculate CoV of expectations
    if i>=nstart
        EU_Mtot_cov=((std(M_tot(1:i,1))/sqrt(i))/EU_Mtot)*100;
        EU_Minert_cov=((std(M_inert(1:i,1))/sqrt(i))/EU_Minert)*100;
    end
    
    %% Increment number of epistemic samples
    i=i+1;
end

end

function [f1,f2]=myEval(fun,u,varargin)
    varargin=varargin{1};
    nVarargs = length(varargin);
    if nVarargs==2
        [f1,f2]=feval(fun,varargin{1},u,varargin{2});
    elseif nVarargs==3
        [f1,f2]=feval(fun,varargin{1},u,varargin{2},varargin{3});
    elseif nVarargs==4
        [f1,f2]=feval(fun,varargin{1},u,varargin{2},varargin{3},varargin{4});
    end
end

