function Env = rlSimulinkEnv(varargin)
% rlSimulinkEnv: Creates RL environment for Simulink models.
%
% ENV = rlSimulinkEnv(MDL,AGENTBLK,OBSERVATIONINFO,ACTIONINFO) Takes a
% Simulink model name MDL, a path to the agent block AGENTBLK, observation
% info OBSERVATIONINFO describing the observation signals, and action info
% ACTIONINFO describing the action signals to create a Simulink environment
% object ENV. OBSERVATIONINFO and ACTIONINFO must be an array of the
% following objects:
%
% 1. <a href="matlab: help rlNumericSpec">rlNumericSpec</a>
% 2. <a href="matlab: help rlFiniteSetSpec">rlFiniteSetSpec</a>
%
% You can also specify a reset function that will be executed before the
% start of each episode during training using the <a href="matlab: help rl.env.SimulinkEnvWithAgent/ResetFcn">ResetFcn</a> property.
%
% Example:
%
%   % setup the arguments
%   mdl = 'rlSimplePendulumModel';
%   agentBlk = [mdl '/RL Agent'];
%   obsInfo = rlNumericSpec([3 1]);
%   obsInfo.Name = 'observations';
%   actInfo = rlNumericSpec([1 1]);
%   actInfo.Name = 'torque';
%
%   % create the environment
%   env = rlSimulinkEnv(mdl,agentBlk,obsInfo,actInfo);
%
%   % include a reset function (randomly initialize theta0 in the model workspace)
%   env.ResetFcn = @(in) setVariable(in,'theta0',randn,'Workspace',mdl);

% Copyright 2018 The MathWorks Inc.

narginchk(4,4);
% catch any errors thrown by the construction of the Simulink environment
% and package them as if they were thrown here
try
    Env = rl.env.SimulinkEnvWithAgent(varargin{:});
catch ex
    throwAsCaller(ex);
end
end