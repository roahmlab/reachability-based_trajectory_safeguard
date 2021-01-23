function Rep = rlStochasticActorRepresentation(Model, ObservationInfo, ActionInfo, varargin)
% RLSTOCHASTICACTORREPRESENTATION construct a stochastic Gaussian actor 
% (policy) representation for reinforcement learning agents.
%
% * Deep Neural Network Representation
%
%       rep = RLSTOCHASTICACTORREPRESENTATION(NET,OINFO,AINFO,'Observation',ONAMES)
%       creates a stochastic gaussian actor representation with default options 
%       using a deep neural network, NET, from Deep Learning Toolbox. 
%       Stochastic actor is typically used in policy gradients agents (e.g.
%       PG, AC, PPO, SAC).
%       Specify the network input layer names, ONAMES, associated with each
%       observation specification as a cell array. The names in ONAMES must
%       be in the same order as the observation specifications, OINFO.
%       Specify the network output layer names, ANAMES, associated with 
%       each action specification as a cell array. The names in ANAMES must
%       be in the same order as the action specifications, AINFO.
%
%       E.g. 2-channel observation, discrete stochastic actor. pi(a1|o) is
%       the probability the actor takes action 1 given observation o.
%                      +-------+
%         o1 --------->|       |          pi(a1|[o1,o2])
%                      |   A   |--------> pi(a2|[o1,o2])
%         o2 --------->|       |          pi(a3|[o1,o2])
%                      +-------+
%       Example of a discrete policy deep neural network representation
%
%       obsInfo = rlNumericSpec([4 1]);
%       actInfo = rlFiniteSetSpec([-10 10]);
%
%       net = [
%           imageInputLayer([4 1 1], 'Normalization', 'none', 'Name', 'state')
%           fullyConnectedLayer(2, 'Name', 'fc')
%           softmaxLayer('Name','actionProb')];
%       actor = rlStochasticActorRepresentation(net, obsInfo, actInfo, 'Observation','state');
%
%       E.g. 2-channel observation, 2-action continuous Gaussian stochastic 
%       actor. REP must have single channel output, which contains the mean 
%       and standard deviation of each continuous action respectively.
%                      +-------+
%         o1 --------->|       |
%                      |   A   |--------> [mu1;mu2;var1;var2]
%         o2 --------->|       |
%                      +-------+
%       Example of a continuous gaussian policy deep neural network representation
%
%       obsInfo = rlNumericSpec([4 1]);
%       actInfo = rlNumericSpec([1 1],'LowerLimit',-10,'UpperLimit',10);
%
%       net = [
%           imageInputLayer([4 1 1], 'Normalization', 'none', 'Name', 'state')
%           fullyConnectedLayer(24, 'Name', 'CommonFC')
%           reluLayer('Name', 'CommonRelu')];
%       meanPath = [
%           fullyConnectedLayer(1, 'Name', 'meanFC')
%           tanhLayer('Name','tanh')
%           scalingLayer('Name','mean','Scale',actInfo.UpperLimit)];
%       stdPath = [
%           fullyConnectedLayer(1, 'Name', 'stdFC')
%           softplusLayer('Name', 'std')
%           ];
%       concatPath = concatenationLayer(3,2,'Name','actionProbParam');
%       net = layerGraph(net);
%       net = addLayers(net, meanPath);
%       net = addLayers(net, stdPath);
%       net = addLayers(net, concatPath);
%       net = connectLayers(net,'CommonRelu','meanFC/in');
%       net = connectLayers(net,'CommonRelu','stdFC/in');
%       net = connectLayers(net,'mean','actionProbParam/in1');
%       net = connectLayers(net,'std','actionProbParam/in2');
%       actor = rlStochasticActorRepresentation(net, obsInfo, actInfo, 'Observation','state');
%
% * Basis Representation
%
%       Policy = RLSTOCHASTICACTORREPRESENTATION({BASISFCN,W0},OINFO,AINFO)
%       creates a basis representation with default options using
%       BASISFCN to define the basis for a policy.
%       B is the column vector returned from B = BASISFCN(obs1,...,obsN),
%       where obs1 to obsN are defined by OINFO.
%       BASISFCN is a function handle.
%       W0 is the initial parameters, a column vector of size [B x (2*NumActions)].
%       - For discrete action space, the model returns action probability 
%         f = softmax(W'*B)
%       - The basis model does not support continuous stochastic action.
%
%   For previous syntax, you can specify nondefault options using an 
%   rlRepresentationOptions object, OPTIONS.
% 
%       REP = RLSTOCHASTICACTORREPRESENTATION(...,OPTIONS)

% Copyright 2019 The MathWorks, Inc.

narginchk(3,6)

% validate data spec
rl.util.validateDataSpecArg(ObservationInfo, mfilename, 'ObservationInfo', 2);
rl.util.validateDataSpecArg(ActionInfo, mfilename, 'ActionInfo', 3);

% parse options
[Options, ObservationNames, ActionNames] = rl.util.parseRepresentationInput(varargin{:});

% neural network specific validations
if ismember(class(Model),{'nnet.cnn.layer.Layer','nnet.cnn.LayerGraph','SeriesNetwork','DAGNetwork','dlnetwork'})
    % neural network model always requires observation names
    if isempty(ObservationNames)
        error(message('rl:agent:errObservationNames'));
    end
    
    % number of obsInfo must be equal to number of obs names
    if numel(ObservationInfo) ~= numel(ObservationNames)
        error(message('rl:agent:errNumObsNamesNeqNumObsInfo'));
    end
end

% basis model specific
if isa(Model,'cell')
    % stochastic actor rep input is observation
    InputSize = {ObservationInfo.Dimension};
    if rl.util.isaSpecType(ActionInfo,'discrete')
        % number of outputs must be equal to number of discrete actions
        NumDiscreteAction = prod(getNumberOfElements(ActionInfo));
        OutputSize = {[NumDiscreteAction 1]};
    else
        % number of output must be 2 * number of continuous actions
        NumContAction = sum(prod(ActionInfo.Dimension),'all');
        OutputSize = {[NumContAction*2 1]};
        
        % NOTE: not support continuous stochastic for basis model
        error(message('rl:agent:errBasisFcnModelContStochasticActorNotSupport'))
    end
else
    InputSize = [];
    OutputSize = [];
end

% create internal model
Model = rl.util.createInternalModelFactory(Model, Options, ObservationNames, ActionNames, InputSize, OutputSize);

Rep = rl.representation.rlStochasticActorRepresentation(...
        Model, ObservationInfo, ActionInfo, Options);