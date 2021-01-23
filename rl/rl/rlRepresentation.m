function rep = rlRepresentation(model,varargin)
% rlRepresentation: Policy and value function representation for reinforcement learning agents
%
% * Deep Neural Network Representation
%
%       rep = rlRepresentation(NET,OINFO,'Observation',ONAMES)
%       creates a representation with default options using NET, which is a
%       deep neural network from Deep Learning Toolbox. Specify the network
%       input layer names, ONAMES, associated with each observation
%       specification as a cell array. The names in ONAMES must be in the
%       same order as the observation specifications, OINFO.
%
%       rep = rlRepresentation(NET,OINFO,AINFO,'Observation',ONAMES,'Action',ANAMES)
%       creates a representation for NET with default options and the
%       specified observation names. Additionally specify the network input
%       layer names, ANAMES, associated with each action specification as a
%       cell array. The names in ANAMES must be the same order as the
%       action specifications, AINFO.
%
% * Table Representation
%
%       rep = rlRepresentation(TABLE)
%       creates a representation with default options using TABLE, which is
%       a table model. To create TABLE, use rlTable.
%
% * Linear Basis Representation
%
%       V = rlRepresentation(BASISFCN,W0,OINFO) 
%       creates a linear basis representation with default options using
%       BASISFCN to define the basis for a value function, such that f =
%       W'*B where B is the column vector returned from B =
%       BASISFCN(obs1,...,obsN), where obs1 to obsN are defined by OINFO.
%       BASISFCN is a function handle.
%
%       Q = rlRepresentation(BASISFCN,W0,{OINFO,AINFO})
%       creates a linear basis representation with default options using
%       BASISFCN to define the basis for a Q function, such that f = W'*B
%       where B is the column vector returned from B =
%       BASISFCN(obs1,...,obsN,act1,...actM), where obs1 to obsN are
%       defined by OINFO and act1 to actM are defined by AINFO. BASISFCN
%       is a function handle.
%
%       Policy = rlRepresentation(BASISFCN,W0,OINFO,AINFO)
%       creates a linear basis representation with default options using
%       BASISFCN to define the basis for a policy, such that f = W'*B
%       where B is the column vector returned from B = BASISFCN(obs1,...,obsN),
%       where obs1 to obsN are defined by OINFO and f corresponds to
%       the actions. BASISFCN is a function handle.
%
%   For all of the previous syntaxes, you can specify nondefault options
%   using an rlRepresentationOptions object, OPTIONS.
% 
%       REP = rlLinearBasisRepresentation(...,OPTIONS)

% Copyright 2018-2019 The MathWorks, Inc.

try
    switch class(model)
        case {'nnet.cnn.layer.Layer','nnet.cnn.LayerGraph','SeriesNetwork','DAGNetwork','dlnetwork'}
            model = rl.internal.dataTransformation.networkToLayerGraph(model);
            rep = iNETRedirect(model,varargin{:});
        case {'rlTable'}
            rep = iTableRedirect(model,varargin{:});
        case {'function_handle'}
            rep = iBasisRedirect(model,varargin{:});
        otherwise
            error(message('rl:agent:errUnknownRepresentation'));
    end
catch ex
    error(message('rl:agent:errObsoleteRepCannotConvertToNewRep'));
end
end

function rep = iTableRedirect(tbl,varargin)
% rep = rlRepresentation(TABLE)
% rep = rlRepresentation(TABLE,OPTS)

[oinfo,ainfo] = getInfo(tbl);
if isempty(varargin)
    opts = rlRepresentationOptions;
else
    opts = varargin{1};
end

if isempty(ainfo)
    % V(o)
    rep = rlValueRepresentation(tbl,oinfo,opts);
else
    % Q(o,a)
    rep = rlQValueRepresentation(tbl,oinfo,ainfo,opts);
end
end

function rep = iBasisRedirect(basisfcn,varargin)
% V = rlRepresentation(BASISFCN,W0,OINFO)
% Q = rlRepresentation(BASISFCN,W0,{OINFO,AINFO})
% A = rlRepresentation(BASISFCN,W0,OINFO,AINFO)

% see if any options are provided
opts = {};
if isa(varargin{end},'rl.option.rlRepresentationOptions')
    opts{1} = varargin{end};
    varargin(end) = [];
end
w0 = varargin{1};

if numel(varargin) > 2
    % A function
    oinfo = varargin{2};
    ainfo = varargin{3};
    % only support discrete stochastic or continuous deterministic actor 
    % in previous releases
    if isa(ainfo,'rl.util.rlNumericSpec')
        rep = rlDeterministicActorRepresentation({basisfcn,w0},oinfo,ainfo,opts);
    else
        rep = rlStochasticActorRepresentation({basisfcn,w0},oinfo,ainfo,opts);
    end
else
    if isa(varargin{end},'cell')
        % Q function
        % only support single action channel in previous releases
        ainfo = varargin{2}{end};
        oinfo = varargin{2}{1:end-1};
        rep = rlQValueRepresentation({basisfcn,w0},oinfo,ainfo,opts);
    else
        % V function
        oinfo = varargin{2};
        rep = rlValueRepresentation({basisfcn,w0},oinfo,opts);
    end
end

end

function rep = iNETRedirect(net,varargin)
% rep = rlRepresentation(NET,OINFO,'Observation',ONAMES)
% rep = rlRepresentation(NET,OINFO,AINFO,'Observation',ONAMES,'Action',ANAMES)
% REVISIT: enable warning in R2020b

% see if any options are provided
opts = {};
if isa(varargin{end},'rl.option.rlRepresentationOptions')
    opts{1} = varargin{end};
    varargin(end) = [];
end

% make sure the net is a layer graph
net = rl.internal.dataTransformation.networkToLayerGraph(net);

if numel(varargin) > 3
    % Q or A function
    anames = varargin{6};
    if all(ismember(anames,net.InputNames))
        % Q function
        % warning(message('rl:agent:warnObsoleteRepConvertToNewRep','rlQValueRepresentation','representation has action inputs'))
        rep = rlQValueRepresentation(net,varargin{:},opts{:});
    else
        % A function
        ainfo = varargin{2};
        if isa(ainfo,'rl.util.rlFiniteSetSpec')
            % discrete
            % warning(message('rl:agent:warnObsoleteRepConvertToNewRep','rlStochasticActorRepresentation','discrete action info is provided'))
            if numel(varargin) > 4
                % detect and remove action names
                ix = find(cellfun(@(x) (ischar(x) || isstring(x)) && strcmpi(x,'Action'),varargin),1);
                varargin(ix:ix+1) = [];
            end
            rep = rlStochasticActorRepresentation(net,varargin{:},opts{:});
        else
            % continuous
            % warning(message('rl:agent:warnObsoleteRepConvertToNewRep','rlDeterministicActorRepresentation','continuous action info is provided'))
            rep = rlDeterministicActorRepresentation(net,varargin{:},opts{:});
        end
    end
else
    % V function
    % rep = RLVALUEREPRESENTATION(NET,OINFO,'Observation',ONAMES)
    % warning(message('rl:agent:warnObsoleteRepConvertToNewRep','rlValueRepresentation','action info is not provided'))
    rep = rlValueRepresentation(net,varargin{:},opts{:});
end
end