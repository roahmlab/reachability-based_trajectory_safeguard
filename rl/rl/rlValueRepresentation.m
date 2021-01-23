function Rep = rlValueRepresentation(Model, ObservationInfo, varargin)
% RLVALUEREPRESENTATION construct state value function representation for
% reinforcement learning agents.
%
% * Deep Neural Network Representation
%
%       rep = RLVALUEREPRESENTATION(NET,OINFO,'Observation',ONAMES)
%       creates a state value representation with default options using a
%       deep neural network, NET, from Deep Learning Toolbox. Specify the
%       network input layer names, ONAMES, associated with each observation
%       specification as a cell array. The names in ONAMES must be in the
%       same order as the observation specifications, OINFO.
%       E.g. 2-channel observation V representation.
%                      +-------+
%         o1 --------->|       |
%                      |   V   |--------> V([o1,o2])
%         o2 --------->|       |
%                      +-------+
%
% * Table Representation
%
%       rep = RLVALUEREPRESENTATION(TABLE,OINFO)
%       creates a state value representation with default options using 
%       TABLE, which is a table model. To create TABLE, use rlTable.
%
% * Basis Representation
%
%       V = RLVALUEREPRESENTATION({BASISFCN,W0},OINFO) 
%       creates a basis representation with default options using
%       BASISFCN to define the basis for a value function, such that f =
%       W'*B where B is the column vector returned from B =
%       BASISFCN(obs1,...,obsN), where obs1 to obsN are defined by OINFO.
%       BASISFCN is a function handle.
%       W0 is the initial parameters, a column vector of size [B x 1].
%
%   For previous syntax, you can specify nondefault options using an
%   rlRepresentationOptions object, OPTIONS.
%
%       REP = RLVALUEREPRESENTATION(...,OPTIONS)

% Copyright 2019 The MathWorks, Inc.

narginchk(2,5)

% validate data spec
rl.util.validateDataSpecArg(ObservationInfo, mfilename, 'ObservationInfo', 2);

% parse options
[Options, ObservationNames, ActionNames] = rl.util.parseRepresentationInput(varargin{:});

% neural network specific validations
if ismember(class(Model),{'nnet.cnn.layer.Layer','nnet.cnn.LayerGraph','SeriesNetwork','DAGNetwork','dlnetwork'})
    % neural network model always requires observation names
    if isempty(ObservationNames)
        error(message('rl:agent:errObservationNames'));
    end
    % number of obsInfo must be equal to number of obsNames
    if numel(ObservationInfo) ~= numel(ObservationNames)
        error(message('rl:agent:errNumObsNamesNeqNumObsInfo'));
    end
end

% basis model specific
if isa(Model,'cell')
    % value rep input is observation
    InputSize = {ObservationInfo.Dimension};
    % value rep output is a scalar
    OutputSize = {[1 1]};
else
    InputSize = [];
    OutputSize = [];
end

% create internal model
Model = rl.util.createInternalModelFactory(Model, Options, ObservationNames, ActionNames, InputSize, OutputSize);

% construct representation
Rep = rl.representation.rlValueRepresentation(...
    Model, ObservationInfo, Options);
end
