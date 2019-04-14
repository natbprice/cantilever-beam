function [Pf,B,MPP,n] = Reliability_Analysis(x,E,opts,model_fidelity)
% Perform reliability analysis using CODES Toolbox
%
% Inputs:
% x  = matrix of design variables
% E = structure with error model
% opts = structure with inputs and options
% model_fidelity = argument specifying model fidelity (optional)
%       low = low-fidelity model 
%       low_upd = low-fidelity model plus Kriging model
%       high = high-fidelity model
%
% Outputs:
% Pf  = probability of failure
% B = reliability index
% MPP = MPP location
% n = safety margin (if requested other outputs set to NaN)
%==========================================================================

%% Unpack variables
dim=opts.reliability.dim;
method=opts.reliability.method;

%% Optional model fidelity specification
if nargin==4
    switch lower(model_fidelity)
        case {'low'}
            e_i=0;
            myFun='g_H_of_E';
            varIn={x,e_i,opts};
        case {'low_upd'}
            myFun='g_H';
            varIn={x,E,opts};
        case {'high'}
            myFun='g_H_true';
            varIn={x,opts};
        otherwise
            error('Unrecognized model fidelity')
    end
else
    myFun='g_H';
    varIn={x,E,opts};
end
g=@(u) myEval(myFun,u,varIn);

%% Calculate reliability or safety margin
if nargout<=3
    Tinv=@(u) myTinv(u,opts);
    
    % Call CODES toolbox to peform reliability calculation
    switch lower(method)
        case {'form'}
            res=CODES.reliability.form(g,dim,'Tinv',Tinv,'display',opts.reliability.form.display,'eps',...
                opts.reliability.form.eps,'solver',opts.reliability.form.solver,'vectorial',opts.reliability.vectorial);
            MPP=res.MPP;
        case {'sorm'}
            res=CODES.reliability.sorm(g,dim,'Tinv',Tinv);
            MPP=res.MPP;
        case {'mcs'}
            res=CODES.reliability.mc(g,dim,'Tinv',Tinv,'CoV',opts.reliability.mcs.cov,'verbose',opts.reliability.mcs.display,...
            'vectorial',opts.reliability.vectorial,'limit',opts.reliability.mcs.nmax,'n',opts.reliability.mcs.nstart);
            MPP=NaN;
        otherwise
            error('Unrecognized reliability method')
    end
    
    Pf=res.Pf;
    B=res.beta;
else
    Pf=NaN; B=NaN; MPP=NaN;
    n=g(opts.design.udet);
end

end

function [g1]=myEval(fun,u,varargin)
    varargin=varargin{1};
    nVarargs = length(varargin);
    if strcmp(fun,'g_H')
        [g1]=feval(fun,varargin{1},u,varargin{2},varargin{3});
    else
        if nVarargs==2
            [g1]=feval(fun,varargin{1},u,varargin{2});
        elseif nVarargs==3
            [g1]=feval(fun,varargin{1},u,varargin{2},varargin{3});
        elseif nVarargs==4
            [g1]=feval(fun,varargin{1},u,varargin{2},varargin{3},varargin{4});
        end
    end
end

