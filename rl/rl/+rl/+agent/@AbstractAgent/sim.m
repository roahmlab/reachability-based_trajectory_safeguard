function varargout = sim(this,env,varargin)
% EXPERIENCES = SIM(ENVIRONMENT,AGENT,OPTIONS) 
% EXPERIENCES = SIM(AGENT,ENVIRONMENT,OPTIONS)
%
% Simulate the ENVIRONMENT against an AGENT with simulation options OPTIONS
% created from rlSimulationOptions.  
%
% An EXPERIENCES structure is returned with the number of elements equal to
% the number of simulations run. The EXPERIENCES structure has the
% following fields:
%
%   Observation         : Structure with fields corresponding to the
%                         observation info of the environment. Each field
%                         is a timeseries of the respective observation
%                         signal. Note, the Observation field is the unit
%                         delay of NextObservation field.
%   Action              : Structure with fields corresponding to the
%                         action info of the environment. Each field
%                         is a timeseries of the respective action
%                         signal.
%   Reward              : timeseries of the scalar reward signal.
%   NextObservation     : See the Observation field.
%   IsDone              : timeseries of the scalar IsDone signal,
%                         signifying the end of episode.
%
% Example: Simulate the Double Integrator environment against an agent for
% up to 1000 steps.
%   
%   env = rlPredefinedEnv('DoubleIntegrator-Discrete'); 
%   opt = rlSimulationOptions('MaxSteps',1000); 
%   experiences = sim(env,agent,opt);
%
% See also: rlSimulationOptions, timeseries

% Copyright 2018 The MathWorks, Inc.

% delegate simulation to environment
[varargout{1:nargout}] = sim(env,this,varargin{:});

end