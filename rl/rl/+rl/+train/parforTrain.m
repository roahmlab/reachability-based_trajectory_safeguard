function parforTrain(trainMgr)
% PARFORTRAIN Uses the training manager to train an agent in parallel using
% parfor. Learning will only happend after the simulations are finished on
% each worker

% Revised: 11-13-2018
% Copyright 2018 The MathWorks, Inc.


agent = trainMgr.Agent;
env = trainMgr.Environment;
trainingOptions = trainMgr.TrainingOptions;
popts = trainingOptions.ParallelizationOptions;

% setup the pool if not already setup
pool = gcp();
numWorkers = pool.NumWorkers;

simOptions = getSimulationOptions(trainingOptions);
% turn off parallel as the simulations will be distributed here with parfor
simOptions.UseParallel = false;
% distribute the env across workers
setupForSim(env,simOptions);
clnup1 = onCleanup(@()cleanupForSim(env));

% set the action msg on the episode mgr
setActionMessage(trainMgr,getString(message('rl:general:TrainingManagerRunningTasksOnWorkers')));
clnupActMsg = onCleanup(@() setActionMessage(trainMgr,getString(message('rl:general:TrainingManagerCleaningUpWorkers'))));

% create worker variables to avoid copies at the start of each parfor
ppoolenv      = parallel.pool.Constant(env  );
ppoolagent    = parallel.pool.Constant(agent);
ppoolmaxsteps = parallel.pool.Constant(trainingOptions.MaxStepsPerEpisode);
ppoolstoperr  = parallel.pool.Constant(trainingOptions.StopOnError);

% the host agent is setup to learn
setStepMode(agent,"learn");
if ~strcmpi(popts.DataToSendFromWorkers,'Gradients')
    % set the worker agents to use exploration and attach loggers to the policy
    parfor i = 1:numWorkers
        % run in exploration mode and attach a logger
        setStepMode(ppoolagent.Value,"sim-with-exploration");
        attachLogger(ppoolagent.Value,ppoolmaxsteps.Value);
    end
    
    setActionMessage(trainMgr,...
        getString(message('rl:general:TrainingManagerRunningTasksOnWorkers')));
    
    % keep training until there are no more remaining sims
    remainingSims = trainingOptions.MaxEpisodes;

    episodeCount = 0;
    while remainingSims > 0 && ~agent.TerminateSimulation
        % figure out how many sims to run (up to the number of available
        % workers)
        activeSims = min(remainingSims,numWorkers);
        workerExperiences = cell(1,activeSims);
        workerEpisodeInfo = cell(1,activeSims);
        simInfos = cell(1,activeSims);
        % get the updated tunable params
        params = getLearnableParameters(agent);
        % run the simulations on the workers
        parfor i = 1:activeSims
            % set the updated tunable parameters
            setLearnableParameters(ppoolagent.Value,params);
            % sim
            [workerExperiences{i},simInfos{i}] = simWithPolicy(ppoolenv.Value,ppoolagent.Value,rlSimulationOptions(...
                'MaxSteps',ppoolmaxsteps.Value,...
                'StopOnError',ppoolstoperr.Value,...
                'NumSimulations',1));
%             % check for term condition
%             workerExperiences{i}{1}{end} = checkForEarlyTermination(ppoolagent.Value,workerExperiences{i}{1}{end});
            workerExperiences{i}{1} = preprocessExperience(ppoolagent.Value, workerExperiences{i}{1});
            % get the episode info from the agent
            workerEpisodeInfo{i} = getEpisodeInfo(ppoolagent.Value);
        end
        % learn from the gathered experiences
        for i = 1:activeSims
            if agent.TerminateSimulation
                break;
            end
            % tell the world data has been received
            notifyDataReceivedFromWorker(trainMgr,rl.comm.TopicMessage("",workerExperiences{i}{1},i));
            % update the episode count
            episodeCount = episodeCount + 1;
            % lean from experiences
            learnFromExperiences(agent,workerExperiences{i}{1});
            % update the training manager
            info.EpisodeInfo = workerEpisodeInfo{i};
            info.SimulationInfo = simInfos{i};
            info.WorkerID = i;
            info.EpisodeCount = episodeCount;
            stopTraining = update(trainMgr,info);
            % terminate the simulation
            terminateSimulation(agent,stopTraining);
            % decrement the remaining sims
            remainingSims = remainingSims - 1;
        end
    end
else
    % set the worker agents to use exploration and attach loggers to the policy
    parfor i = 1:numWorkers
        % run in exploration mode, append experience and attach a logger
        setStepMode(ppoolagent.Value,"sim-with-exploration-and-append");
        attachLogger(ppoolagent.Value,ppoolmaxsteps.Value);
    end
    
    setActionMessage(trainMgr,...
        getString(message('rl:general:TrainingManagerRunningTasksOnWorkers')));
    
    % keep training until there are no more remaining sims
    remainingSims = trainingOptions.MaxEpisodes;
%     simInfos = cell(1,remainingSims);
    episodeCount = 0;
    while remainingSims > 0 && ~agent.TerminateSimulation
        % figure out how many sims to run (up to the number of available
        % workers)
        activeSims = min(remainingSims,numWorkers);
        workerGradients = cell(1,activeSims);
        workerEpisodeInfo = cell(1,activeSims);
        simInfos = cell(1,activeSims);
        % get the updated tunable params
        params = getLearnableParameters(agent);
        % run the simulations on the workers
        parfor i = 1:activeSims
            % set the updated tunable parameters
            setLearnableParameters(ppoolagent.Value,params);
            % sim
            [~,simInfos{i}] = simWithPolicy(ppoolenv.Value,ppoolagent.Value,rlSimulationOptions(...
                'MaxSteps',ppoolmaxsteps.Value,...
                'StopOnError',ppoolstoperr.Value,...
                'NumSimulations',1));
            % get the episode info from the agent
            workerEpisodeInfo{i} = getEpisodeInfo(ppoolagent.Value);
            % accumulate gradients from learned episode
            workerGradients{i} = accumulateGradient(ppoolagent.Value,workerEpisodeInfo{i}.StepsTaken);
        end
        % learn from the gathered gradients
        for i = 1:activeSims
            if agent.TerminateSimulation
                break;
            end
            % tell the world data has been received
            notifyDataReceivedFromWorker(trainMgr,rl.comm.TopicMessage("",workerGradients{i},i));
            % update the episode count
            episodeCount = episodeCount + 1;
            % apply gradients from workers
            applyGradient(agent,workerGradients{i});
            % update the training manager
            info.EpisodeInfo = workerEpisodeInfo{i};
            info.SimulationInfo = simInfos{i};
            info.WorkerID = i;
            info.EpisodeCount = episodeCount;
            stopTraining = update(trainMgr,info);
            % terminate the simulation
            terminateSimulation(agent,stopTraining);
            % decrement the remaining sims
            remainingSims = remainingSims - 1;
        end
    end
end
% % attach the sim info
% trainMgr.SimulationInfo = vertcat(simInfos{:});

