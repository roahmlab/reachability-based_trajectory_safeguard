function Options = rlTrainingOptions(varargin)
% RLTRAININGOPTIONS: Creates training options for RL.
%
%   OPT = RLTRAININGOPTIONS returns the default options for RL agents. 
%
%   OPT = RLTRAININGOPTIONS('Option1',Value1,'Option2',Value2,...) uses name/value
%   pairs to override the default values for 'Option1','Option2',...
%
%   Supported options are:
%
%   MaxEpisodes                 Maximum number of episodes to train the agent
%   MaxStepsPerEpisode          Maximum number of steps to run per episode
%   ScoreAveragingWindowLength  Window length for averaging scores and rewards
%   StopTrainingCriteria        Stop criteria for agent training
%                               - AverageSteps   : Running average of number of steps per episode
%                               - AverageReward  : Running average of reward per episode
%                               - EpisodeReward  : Reward for current episode
%                               - GlobalStepCount: Total times the agent was invoked
%                               - EpisodeCount   : Total number of episodes the agent has trained for
%   StopTrainingValue           Value of training criteria to stop
%   SaveAgentCriteria           Criteria to save agent while training
%                               - None           : Do not save any agent in training
%                               - AverageSteps   : Running average of number of steps per episode
%                               - AverageReward  : Running average of reward per episode
%                               - EpisodeReward  : Reward for current episode
%                               - GlobalStepCount: Total times the agent was invoked
%                               - EpisodeCount   : Total number of episodes the agent has trained for
%   SaveAgentValue              Value of criteria to save
%   SaveAgentDirectory          Directory name for saving agents
%   UseParallel                 Use parallel training schemes defined by ParallelizationOptions. Otherwise the agent is trained on a single process (false (default), true)
%   ParallelizationOptions      Parallelization options to control parallel training schemes
%                               - Mode                           : Parallelization scheme for training
%                                   - sync (default)             : Use the parpool to run synchronous training on the available workers
%                                   - async                      : Use the parpool to run asynchronous training on the available workers
%                               - DataToSendFromWorkers          : Data that will be sent from the workers to the client ("experiences" (default), "gradients").
%                               - StepsUntilDataIsSent           : The number of steps the worker will take before sending data back to the client. If -1 (default), the worker will send data at the end of an episode.
%                               - WorkerRandomSeeds              : Initialize the random seed on the worker using <a href="matlab: help rng">rng</a>. If -1 (default) the random seed on each worker will be initialized with the worker ID. Use a vector to initialize a seed for each worker.
%                               - TransferBaseWorkspaceVariables : Transfer base workspace variables to each worker before training starts ("on" (default), "off"). 
%                               - AttachedFiles                  : Files to attach to the parpool
%                               - SetupFcn                       : Setup function handle to run on each worker before training starts.
%                               - CleanupFcn                     : Cleanup function handle to run on each worker after training ends.
%   StopOnError                 Stop remaining simulations if any simulation errors occurs. ("on" (default), "off")
%   Verbose                     Display training progress on the command line (false (default), true)
%   Plots                       Display training progress with the Episode Manager ("training-progress" (default), "none")
%
% See also: <a href="matlab: help rl.agent.AbstractAgent\train">train</a>, parpool

% Copyright 2017-2018 The MathWorks Inc.

Options = rl.option.rlTrainingOptions(varargin{:});

end