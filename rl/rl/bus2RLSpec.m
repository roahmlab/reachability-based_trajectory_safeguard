function specs = bus2RLSpec(varargin)
% SPECS = BUS2RLSPEC(BUSNAME)
%
% Create a set of RL data specifications from a Simulink bus. Each leaf
% element in the bus will correspond to an RL data specification element.
% The bus must be defined in the base workspace.
%
% SPECS = BUS2RLSPEC(BUSNAME,NAME1,VALUE1,...)
% 
% Provide additional name-value pairs to control how the RL data
% specifications are created. The supported name-value pairs are:
%
%   Model               : Provide the model name if the bus is defined in
%                         the model global workspace (e.g. Data
%                         Dictionary) rather than the base workspace.
%   BusElementNames     : The element names to generate RL data
%                         specifications for. The provided elements must be
%                         leaf elements of the bus. By default, RL data
%                         specifications will be created for all leaf bus
%                         elements. BusElementNames should not be specified
%                         when creating specifications for action signals
%                         since the RL Agent block must output the full bus
%                         signal.
%   DiscreteElements    : Specify a finite set element for a specific bus
%                         element via name-value pair cell array. The
%                         specified discrete element set must be castable
%                         to the element data type.
%
% Example 1: Create RL data specifications for an observation set from a
% bus in the base workspace.
%
%   mdl = 'SimplePendulumModel_bus'
%
%   % create the observation bus
%   obsBus = Simulink.Bus();
%   elo(1) = Simulink.BusElement;
%   elo(1).Name = 'sin_theta';
%   elo(2) = Simulink.BusElement;
%   elo(2).Name = 'cos_theta';
%   elo(3) = Simulink.BusElement;
%   elo(3).Name = 'dtheta';
%   obsBus.Elements = elo;
%
%   % create the data spec
%   observations = bus2RLSpec('obsBus');
%
% Example 2: Create RL data specifications for an observation set from a
% bus defined in a data dictionary.
%
%   observations = bus2RLSpec('obsBus','Model',mdl);
%
% Example 3: Create RL data specifications for the first two elements of
% the bus.
%
%   observations = bus2RLSpec('obsBus','BusElementNames',["sin_theta","cos_theta"]);
%
% Example 4: Create RL data specifications where the first two elements are
% n-dimensional signals and the third element is a finite set with values =
% [-1 0 1]
%
%   observations = bus2RLSpec('obsBus','DiscreteElements',{'dtheta',-1:1});
%
% See also: rlSimulinkEnv, rlNumericSpec, rlFiniteSetSpec

% Revised: 12-3-2018
% Copyright 2018 The MathWorks, Inc.

try
    specs = rl.util.bus2RLData(varargin{:});
catch ex
    throwAsCaller(ex);
end


