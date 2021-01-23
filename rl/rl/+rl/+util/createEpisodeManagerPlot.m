function episodeManager = createEpisodeManagerPlot(trainingResult)
% Create Episode manager plot with trainingResult struct. 
% 
% trainingResult can be either result from saved agent or from stopped
% training.
%
% trainingResult from saved agent has the following structure:
%            savedAgentInfo = struct(...
%                 'TrainingOpts',                 TrainingOptions,...
%                 'HasCritic',                    HasCritic,...
%                 'HardwareResource',             Devices,...
%                 'LearningRate',                 LearnRates,...
%                 'TraningStartTime',             TrainingStartTime,...
%                 'ElapsedTime',                  elapsedTime);
%             savedAgentResultStruct = struct(...
%                 'TrainingStats',                stats,...
%                 'Information',                  savedAgentInfo);           
%
% trainingResult from stopped training has the following structure:
%             trainingInfo = struct(...
%                 'TrainingOpts',                 TrainingOptions,...
%                 'HasCritic',                    HasCritic,...
%                 'HardwareResource',             Devices,...
%                 'LearningRate',                 LearnRates,...
%                 'TraningStartTime',             TrainingStartTime,...
%                 'ElapsedTime',                  TrainingElapsedTime,...
%                 'StopTrainingCriteria',         trainingStoppedReason,...
%                 'StopTrainingValue',            trainingStoppedValue...
%                 );
%              trainingResultStruct = struct(...
%                 'EpisodeIndex',                 vector1, ...
%                 'EpisodeReward',                vector1, ...
%                 'EpisodeSteps',                 vector1, ...
%                 'AverageReward',                vector1, ...
%                 'AverageSteps',                 vector1, ...
%                 'TotalAgentSteps',              vector1, ...
%                 'Information',                  trainingInfo);     

% Copyright 2018-2019 The MathWorks, Inc.

%% 
% create episode manager plot.
res = trainingResult;
% traaining result info
info = res.Information;
% create episode manager with data to be updated.
episodeManager = rl.util.EpisodeManager(info.TrainingOpts,info.HardwareResource,info.LearningRate,info.HasCritic);
% update epsidoe manager with date from trainingResult.
if  isfield(info,'StopTrainingCriteria')
    postTraining = true;
else
    postTraining = false;
end
updateWithTrainingResult(episodeManager,res,postTraining);
end