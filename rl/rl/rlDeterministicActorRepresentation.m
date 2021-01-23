function Rep = rlDeterministicActorRepresentation(Model, ObservationInfo, ActionInfo, varargin)
% RLDETERMINISTICACTORREPRESENTATION construct a deterministic actor 
% (policy) representation for reinforcement learning agents.
%
% * Deep Neural Network Representation
%
%       REP = RLDETERMINISTICACTORREPRESENTATION(NET,OINFO,AINFO,'Observation',ONAMES,'Action',ANAMES)
%       creates a deterministic actor representation with default options 
%       using a deep neural network, NET, from Deep Learning Toolbox. 
%       Deterministic actor is used in continuous action space problems.
%       Specify the network input layer names, ONAMES, associated with each
%       observation specification as a cell array. The names in ONAMES must
%       be in the same order as the observation specifications, OINFO.
%       Specify the network output layer names, ANAMES, associated with 
%       each action specification as a cell array. The names in ANAMES must
%       be in the same order as the action specifications, AINFO.
%       E.g. 2-channel observation, 2-channel deterministic action actor.
%                      +-------+
%         o1 --------->|       |--------> a1
%                      |   A   |
%         o2 --------->|       |--------> a2
%                      +-------+
%
% * Basis Representation
%
%       Policy = RLDETERMINISTICACTORREPRESENTATION({BASISFCN,W0},OINFO,AINFO)
%       creates a basis representation with default options using
%       BASISFCN to define the basis for a policy, such that f = W'*B
%       where B is the column vector returned from B = BASISFCN(obs1,...,obsN),
%       where obs1 to obsN are defined by OINFO and f corresponds to
%       the actions. 
%       BASISFCN is a function handle.
%       W0 is the initial parameters, a column vector of size [B x NumActions].
%
%   For previous syntax, you can specify nondefault options using an 
%   rlRepresentationOptions object, OPTIONS.
% 
%       REP = RLDETERMINISTICACTORREPRESENTATION(...,OPTIONS)

% Copyright 2019 The MathWorks, Inc.

narginchk(3,8)

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
    
    % neural network model always requires action names
    if isempty(ActionNames)
        error(message('rl:agent:errActionNames'));
    end
    
    % number of actInfo must be equal to number of act names
    if numel(ActionInfo) ~= numel(ActionNames)
        error(message('rl:agent:errNumActNamesNeqNumActInfo'));
    end
end

% basis model specific
if isa(Model,'cell')
    % deterministic actor rep input is observation
    InputSize = {ObservationInfo.Dimension};
    % deterministic actor rep output is action
    OutputSize = {ActionInfo.Dimension};
else
    InputSize = [];
    OutputSize = [];
end

% create internal model
Model = rl.util.createInternalModelFactory(Model, Options, ObservationNames, ActionNames, InputSize, OutputSize);

Rep = rl.representation.rlDeterministicActorRepresentation(...
        Model, ObservationInfo, ActionInfo, Options);