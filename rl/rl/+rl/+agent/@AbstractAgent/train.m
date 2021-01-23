function TrainingStatistics = train(this,env,varargin)
% TRAININGSTATISTICS = TRAIN(AGENT,ENVIRONMENT,OPTIONS)
%
% Train the AGENT against the ENVIRONMENT with training options OPTIONS
% created from rlTrainingOptions. 
%
% TRAININGSTATISTICS are returned when training is finished detailing the
% history of the training progress. TRAININGSTATISTICS is a structure with
% the following fields:
%
%   EpisodeIndex        : Episode number.
%   EpisodeReward       : Cumulative episode reward.
%   EpisodeSteps        : Steps taken during episode.
%   AverageReward       : Average cumulative reward of episodes in the
%                         averaging window (specified by
%                         ScoreAveragingWindowLength in the training
%                         options).
%   TotalAgentSteps     : Total number of steps taken.
%   EpisodeQ0           : Critic estimate of the long-term reward at the
%                         initial conditions of the environment (only for
%                         agents with a critic e.g. DQN and DDPG).
%
% See also: rlTrainingOptions, <a href="matlab: help rl.env.Abstract\sim">sim</a>

% Copyright 2018 The MathWorks, Inc.

parser = inputParser();
addRequired(parser,'Environment',@(x)isa(x,'rl.env.AbstractEnv'));
addOptional(parser,'TrainingOptions',rlTrainingOptions(),@(x)isa(x,'rl.option.rlTrainingOptions'));
parse(parser,env,varargin{:});

trainingOptions = parser.Results.TrainingOptions;
validateStopTrainingFunction(trainingOptions);
validateSaveAgentFunction(trainingOptions);
trainingOptions = validateAgentTrainingCompatibility(this,trainingOptions);
% set the options
this.TrainingOptions = trainingOptions;

% validate observation and action info of env and agent
envActionInfo = getActionInfo(env);
envObservationInfo = getObservationInfo(env);
if ~isCompatible(this.ActionInfo,envActionInfo)
    error(message('rl:general:errActionInfo'));
end
if ~isCompatible(this.ObservationInfo,envObservationInfo)
    error(message('rl:general:errObservationInfo'));    
end

%% create the training manager
trainMgr = rl.train.TrainingManager(env,this,trainingOptions);
clnup = onCleanup(@() delete(trainMgr));

%% run the training
TrainingStatistics = run(trainMgr);