function Rep = rlQValueRepresentation(Model, ObservationInfo, ActionInfo, varargin)
% RLQVALUEREPRESENTATION construct state-action value function 
% representation for reinforcement learning agents.
%
% * Deep Neural Network Representation
%
%       rep = RLQVALUEREPRESENTATION(NET,OINFO,AINFO,'Observation',ONAMES)
%       creates a multi-output state-action value representation for 
%       discrete action space with default options using a deep neural 
%       network, NET, from Deep Learning Toolbox.
%       NET must have single output layer with number of outputs equal to 
%       number of possible discrete actions.
%       Specify the network input layer names, ONAMES, associated with each 
%       observation specification as a cell array. The names in ONAMES must
%       be in the same order as the observation specifications, OINFO.
%       E.g. 2-channel observation, 3 discrete actions Q representation.
%                      +-------+
%         o1 --------->|       |          Q([o1,o2],a1)
%                      |   Q   |--------> Q([o1,o2],a2)
%         o2 --------->|       |          Q([o1,o2],a3)
%                      +-------+
%
%       rep = RLQVALUEREPRESENTATION(NET,OINFO,AINFO,'Observation',ONAMES,'Action',ANAMES)
%       creates a single-output state-action value representation with 
%       default options for NET. Additionally specify the network input
%       layer names, ANAMES, associated with each action specification as a
%       cell array. The names in ANAMES must be the same order as the
%       action specifications, AINFO.
%       E.g. 2-channel observation, single channel action Q single-output 
%       representation.
%                      +-------+
%         o1 --------->|       |
%         o2 --------->|   Q   |--------> Q([o1,o2],a)
%         a  --------->|       |
%                      +-------+
%
% * Table Representation
%
%       rep = RLQVALUEREPRESENTATION(TABLE,OINFO,AINFO)
%       creates a state value representation with default options using 
%       TABLE, which is a table model. To create TABLE, use rlTable.
%
% * Basis Representation
%
%       Q = RLQVALUEREPRESENTATION({BASISFCN,W0},OINFO,AINFO)
%       creates a basis representation with default options using
%       BASISFCN to define the basis for a Q function, such that f = W'*B
%       where B is the column vector returned from B =
%       BASISFCN(obs1,...,obsN,act1,...actM), where obs1 to obsN are
%       defined by OINFO and act1 to actM are defined by AINFO. 
%       BASISFCN is a function handle.
%       W0 is the initial parameters, a column vector of size [B x 1].
%
%   For all of the previous syntaxes, you can specify nondefault options
%   using an rlRepresentationOptions object, OPTIONS.
% 
%       REP = RLQVALUEREPRESENTATION(...,OPTIONS)

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
    
    % number of actInfo must be equal to number of act names
    if ~isempty(ActionNames) && numel(ActionInfo) ~= numel(ActionNames)
        error(message('rl:agent:errNumActNamesNeqNumActInfo'));
    end
end

% basis model specific validations
if isa(Model,'cell')
    BasisFunction = Model{1};
    if ~isa(BasisFunction,'function_handle')
        error(message('rl:agent:errInvalidInputBasisModel'));
    end
    NumModelInputChannel = nargin(BasisFunction);
    switch NumModelInputChannel
        case numel(ObservationInfo)
            % multi output Q
            InputSize = {ObservationInfo.Dimension};
            if rl.util.isaSpecType(ActionInfo, 'continuous')
                error(message('rl:agent:errQContinuousActionMultiOutput'));
            end
            OutputSize = {[getNumberOfElements(ActionInfo) 1]};
        case numel(ObservationInfo) + numel(ActionInfo)
            % single output Q
            InputSize = {ObservationInfo.Dimension ActionInfo.Dimension};
            OutputSize = {[1 1]};
        otherwise
            error(message('rl:agent:errInvalidNumInputBasisModel'));
    end
else
    InputSize = [];
    OutputSize = [];
end

% create internal model
Model = rl.util.createInternalModelFactory(Model, Options, ObservationNames, ActionNames, InputSize, OutputSize);

% check action names are provided for single output Q (neural network)
if isa(Model,'rl.representation.model.rlLayerModel')
    OutputSize = getSize(Model,'output');
    NumOutput = prod(OutputSize{1});
    if NumOutput == 1
        if isempty(ActionNames)
            error(message('rl:agent:errSingleQLayerActNamesMissing'));
        end
    end
end

% construct representation
Rep = rl.representation.rlQValueRepresentation(Model, ...
    ObservationInfo, ActionInfo, Options);
end
