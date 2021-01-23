function env = rlPredefinedEnv(keyword)
% RLPREDEFINEDENV: Creates predefined RL environment.
%
% ENV = RLPREDEFINEDENV(KEYWORD) Creates the predefined RL environment
% ENV that corresponds to the specified keyword KEYWORD. You can create
% both MATLAB and Simulink predefined environments.
%
% Supported MATLAB environment keywords are:
%
%   'BasicGridWorld'
%   'CartPole-Discrete'
%   'CartPole-Continuous'
%   'DoubleIntegrator-Discrete'
%   'DoubleIntegrator-Continuous'
%   'SimplePendulumWithImage-Discrete'
%   'SimplePendulumWithImage-Continuous'
%   'WaterFallGridWorld-Deterministic'
%   'WaterFallGridWorld-Stochastic'
%
% Supported Simulink environment keywords are:
%
%   'SimplePendulumModel-Discrete'
%   'SimplePendulumModel-Continuous'
%   'CartPoleSimscapeModel-Discrete'
%   'CartPoleSimscapeModel-Continuous'
%
% Examples:
%   
%   % create discrete cart pole MATLAB environment
%   env = rlPredefinedEnv('CartPole-Discrete');
%
%   % create continuous cart pole MATLAB environment
%   env = rlPredefinedEnv('CartPole-Continuous');
%   
%   % create discrete pendulum Simulink environment
%   env = rlPredefinedEnv('SimplePendulumModel-Discrete');
%
%   % create continuous pendulum Simulink environment
%   env = rlPredefinedEnv('SimplePendulumModel-Continuous');

% Copyright 2018 The MathWorks, Inc.

validateattributes(keyword,{'char','string'},{'scalartext'})
keyword = lower(char(keyword));

switch keyword
    case lower('BasicGridWorld')
        env = rl.env.BasicGridWorld;
    case lower('WaterFallGridWorld-Deterministic')
        env = rl.env.WaterFallDeterministicGridWorld;
    case lower('WaterFallGridWorld-Stochastic')
        env = rl.env.WaterFallStochasticGridWorld;
    case lower('CartPole-Discrete')
        env = rl.env.CartPoleDiscreteAction;
    case lower('CartPole-Continuous')
        env = rl.env.CartPoleContinuousAction;  
    case lower('DoubleIntegrator-Discrete')
        env = rl.env.DoubleIntegratorDiscreteAction;
    case lower('DoubleIntegrator-Continuous')
        env = rl.env.DoubleIntegratorContinuousAction;
    case lower('SimplePendulumWithImage-Discrete')
        env = rl.env.SimplePendlumWithImageDiscreteAction;
    case lower('SimplePendulumWithImage-Continuous')
        env = rl.env.SimplePendlumWithImageContinuousAction;
    case lower('SimplePendulumModel-Discrete')
        env = localCreateSimplePendulumEnv('discrete');
    case lower('SimplePendulumModel-Continuous')
        env = localCreateSimplePendulumEnv('continuous');
    case lower('CartPoleSimscapeModel-Discrete')
        env = localCreateCartPoleSimscapeEnv('discrete');
    case lower('CartPoleSimscapeModel-Continuous')
        env = localCreateCartPoleSimscapeEnv('continuous');    
    otherwise
        error(message('rl:env:rlPredefinedEnvInvalidKeyword'));
end

function env = localCreateSimplePendulumEnv(actionType)
% crete env from simple pendulum model
mdl = 'rlSimplePendulumModel';
blk = 'rlSimplePendulumModel/RL Agent';
load_system(mdl);
if strcmp(get_param('rlSimplePendulumModel/create observations','ThetaObservationHandling'),'sincos')
    numObs = 3;
else
    numObs = 2;
end

% create the observation info
observationInfo = rlNumericSpec([numObs 1]);
observationInfo.Name = 'observations';

% create the action info
if strcmp(actionType,'continuous')
    actionInfo = rlNumericSpec([1 1],...
        'LowerLimit',-2,'UpperLimit',2);
else
    actionInfo = rlFiniteSetSpec([-2 0 2]);
end
actionInfo.Name = 'torque';

env = rlSimulinkEnv(mdl,blk,observationInfo,actionInfo);

function env = localCreateCartPoleSimscapeEnv(actionType)
% create env from cart pole simscape model
mdl = 'rlCartPoleSimscapeModel';
blk = 'rlCartPoleSimscapeModel/RL Agent';
load_system(mdl);

% create the observation info
observationInfo = rlNumericSpec([5 1]);
observationInfo.Name = 'observations';

% create the action info
if strcmp(actionType,'continuous')
    actionInfo = rlNumericSpec([1 1],...
        'LowerLimit',-15,'UpperLimit',15);
else
    actionInfo = rlFiniteSetSpec([-15 0 15]);
end
actionInfo.Name = 'force';

env = rlSimulinkEnv(mdl,blk,observationInfo,actionInfo);
