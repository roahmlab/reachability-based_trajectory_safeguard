function Options = rlSimulationOptions(varargin)
% RLSIMULATIONOPTIONS: Creates simulation options.
%
%   OPT = RLSIMULATIONOPTIONS returns the default simulation options for RL
%   environment and agent.
%
%   OPT = RLSIMULATIONOPTIONS('Option1',Value1,'Option2',Value2,...) uses name/value
%   pairs to override the default values for 'Option1','Option2',...
%
%   Supported options are:
%
%   MaxSteps                    Number of steps for simulation
%   NumSimulations              Number of simulations to run
%   StopOnError                 Stop remaining simulations if any simulation errors occur. ("on" (default), "off")
%   UseParallel                 Run simulations in parallel. Otherwise the simulation is run on a single process (false (default), true)
%   ParallelizationOptions      Parallelization options to control parallel simulation
%                               - WorkerRandomSeeds              : Initialize the random seed on the worker using <a href="matlab: help rng">rng</a>. If -1 (default) the random seed on each worker will be initialized with the worker ID. Use a vector to initialize a seed for each worker.
%                               - TransferBaseWorkspaceVariables : Transfer base workspace variables to each worker before simulations start ("on" (default), "off"). 
%                               - AttachedFiles                  : Files to attach to the parpool
%                               - SetupFcn                       : Setup function handle to run on each worker before simulations start.
%                               - CleanupFcn                     : Cleanup function handle to run on each worker after simulations ends.
%
% See also: <a href="matlab: help rl.env.Abstract\sim">sim</a>

% Copyright 2017-2018 The MathWorks Inc.

Options = rl.option.rlSimulationOptions(varargin{:});

end