function seriesTrain(trainMgr)
% SERIESTRAIN
% Perform series training on a local workstation

% Revised: 10-25-2018
% Copyright 2018 The MathWorks, Inc.

agent = trainMgr.Agent;
env = trainMgr.Environment;
trainingOptions = trainMgr.TrainingOptions;

%% initialize
setStepMode(agent,"learn");

%% run the outer episode loop
[~,simInfo] = simWithPolicy(env,agent,rlSimulationOptions(...
        'MaxSteps',trainingOptions.MaxStepsPerEpisode,...
        'NumSimulations',trainingOptions.MaxEpisodes,...
        'StopOnError',trainingOptions.StopOnError));
    
% attach the simulation info to the training manager
trainMgr.SimulationInfo = simInfo;



