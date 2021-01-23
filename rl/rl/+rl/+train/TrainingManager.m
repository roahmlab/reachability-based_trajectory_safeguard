classdef TrainingManager < handle
% TRAININGMANAGER

% Revised: 11-2-2018
% Copyright 2018 The MathWorks, Inc.

    properties
        Agent 
        Environment 
        TrainingOptions 
    end
    properties (Hidden, SetAccess = private)
        
        HasCritic (1,1) logical
        
        % episode count "states"
        EpisodeCount = 0
        TotalEpisodeStepCount = 0
        
        % window "states"
        RewardsInAveragingWindow
        StepsInAveragingWindow
        
        % construct training result struct
        TrainingStartTime
        TrainingElapsedTime
        Watch
        Devices
        LearnRates

		% simulation info
        SimulationInfo = {}
        
        TrainingStats
    end
    properties (Access = private,Transient)
        % listener to determine if stop training has been requested
        StopTrainingListener
        
        % listener for episode finished event
        EpisodeFinishedListener
        
        % listener for data rcv on worker
        DataReceivedOnWorkerListener
        
        % listener for tasks on workers
        TasksRunningOnWorkersListener
        TasksCleanedUpOnWorkersListener
        TasksActionMessageReceived
        
        
        % episode manager
        EpisodeMgr
        
        % pre train settings
        PreTrainAgentSettings = []
    end
    events
        TrainingManuallyStopped
        DataReceivedFromWorker
        TrainingManagerUpdated
    end
    methods
        function this = TrainingManager(env,agent,opt)
            this.Environment = env;
            this.Agent = agent;
            this.TrainingOptions = opt;
            this.HasCritic = hasCritic(this.Agent);
            
            % watch
            this.Watch = nnet.internal.cnn.ui.adapter.Stopwatch();
            reset(this.Watch);
            
            % training start time
            this.TrainingStartTime = iDateTimeAsStringUsingDefaultLocalFormat(datetime('now'));
        end
        function delete(this)
            delete(this.StopTrainingListener);
            delete(this.EpisodeFinishedListener);
            delete(this.DataReceivedOnWorkerListener);
            delete(this.TasksRunningOnWorkersListener);
            delete(this.TasksCleanedUpOnWorkersListener);
            delete(this.TasksActionMessageReceived);
        end 
        function cleanup(this)
            setStepMode(this.Agent,"sim");
            % only call postTrain if preTrain was successfully called,
            % which will return a struct.
            if ~isempty(this.PreTrainAgentSettings)
                postTrain(this.Agent,this.PreTrainAgentSettings);
            end
            this.PreTrainAgentSettings = [];
            terminateSimulation(this.Agent,false);
        end              
        function setActionMessage(this,msg)
            % set an action message for the episode manager
            if isvalid(this)
                episodeMgr = this.EpisodeMgr;
                if ~isempty(episodeMgr)
                    setActionMessage(episodeMgr,msg);
                end
            end
        end
        function msg = getActionMessage(this)
            % get an action message from the episode manager
            msg = '';
            if isvalid(this)
                episodeMgr = this.EpisodeMgr;
                if ~isempty(episodeMgr)
                    msg = getActionMessage(episodeMgr);
                end
            end
        end
        function reqSimulink = requiresSimulink(this)
            % is simulink needed to run the training
            reqSimulink = isa(this.Environment,'rl.env.SimulinkEnvWithAgent');
        end
        function stopTraining = update(this,episodeFinishedInfo)
            % update the manager once an episode finishes
            
            epinfo       = episodeFinishedInfo.EpisodeInfo   ;
            episodeCount = episodeFinishedInfo.EpisodeCount  ;
            workerID     = episodeFinishedInfo.WorkerID      ;
            simInfo      = episodeFinishedInfo.SimulationInfo;
            
            % DEBUG
            % fprintf('finished sim C %u\n',episodeCount);
            
            % attach the episode info
            this.SimulationInfo{episodeCount} = simInfo;
            
            % update "states"
            this.EpisodeCount = episodeCount;
            this.TotalEpisodeStepCount = this.TotalEpisodeStepCount + epinfo.StepsTaken;
            
            % compute info and update displays
            info = computeMetrics(this,epinfo);
            stopTraining = updateDisplaysFromTrainingInfo(this,info);
            
            % terminate the simulation if criteria is met
            terminateSimulation(this.Agent,stopTraining);
            
            % tell the world that the training manager has been updated
            s = struct;
            s.EpisodeInfo = epinfo;
            s.ComputedEpisodeInfo = info;
            s.EpisodeCount = episodeCount;
            s.WorkerID = workerID;
            ed = rl.util.RLEventData(s);
            notify(this,'TrainingManagerUpdated',ed);
            
            % DEBUG
            % fprintf('train mgr update on id = %u\n',workerID);
        end
        function stats = run(this)
            % run the training routine with setup and cleanup
            clnup = onCleanup(@() cleanup(this));
            
            preTrain(this);
            train(this);
            stats = postTrain(this);
        end
        function preTrain(this)
            % run before training occurs
            
            % make sure the agent is reset
            agent = this.Agent;
            agent.MaxSteps = this.TrainingOptions.MaxStepsPerEpisode;
            this.PreTrainAgentSettings = preTrain(agent);
            terminateSimulation(agent,false);
            
            % initialize the window "states"
            numScoresToAverage = this.TrainingOptions.ScoreAveragingWindowLength;
            this.RewardsInAveragingWindow = zeros(numScoresToAverage,1,'single');
            this.StepsInAveragingWindow   = zeros(numScoresToAverage,1,'single');
            
            
            % build the train stats struct
            maxEpisodes = this.TrainingOptions.MaxEpisodes;
            this.SimulationInfo = cell(1,maxEpisodes);
            % build the saved agents directory
            if isempty(dir(this.TrainingOptions.SaveAgentDirectory))
                if ~strcmpi(this.TrainingOptions.SaveAgentCriteria,"none")
                    try
                        mkdir(this.TrainingOptions.SaveAgentDirectory);
                    catch ex
                        me = MException(message('rl:general:TrainingManagerUnableToCreateSaveDir',this.TrainingOptions.SaveAgentDirectory));
                        throw(addCause(me,ex));
                    end
                end
            end
            
            % Initialize training statistics
            TrainingStatistics = struct(...
                'EpisodeIndex',zeros(maxEpisodes,1), ...
                'EpisodeReward',zeros(maxEpisodes,1), ...
                'EpisodeSteps',zeros(maxEpisodes,1), ...
                'AverageReward',zeros(maxEpisodes,1), ...
                'AverageSteps',zeros(maxEpisodes,1),...
                'TotalAgentSteps',zeros(maxEpisodes,1),...
                'Information',[]);
            if this.HasCritic
                TrainingStatistics.EpisodeQ0 = zeros(maxEpisodes,1);
            end
            this.TrainingStats = TrainingStatistics;
            % create the episode manager
            initializeEpisodeManager(this);
        end
        function stats = postTrain(this)
            % return the training stats post train     
            maxEpisodes = this.TrainingOptions.MaxEpisodes;
            episodeIndex = this.EpisodeCount;
            this.TrainingElapsedTime = char(getDurationSinceReset(this.Watch));
            % determine the stopping reason and value
            if episodeIndex == maxEpisodes
                % TODO MSG
                trainingStoppedReason = 'Maximum number of episodes';
                trainingStoppedValue = episodeIndex;
            else
                trainingStoppedReason = this.TrainingOptions.StopTrainingCriteria;
                trainingStoppedValue =  this.TrainingOptions.StopTrainingValue;
            end
            
            % tell the episode manager that training has stopped
            episodeMgr = this.EpisodeMgr;
            if ~isempty(episodeMgr) && isvalid(episodeMgr)
                setElapsedTime(this.EpisodeMgr,this.TrainingElapsedTime);
                stopTraining(episodeMgr,trainingStoppedReason,trainingStoppedValue);
            end
            
            % Clean up unused training statistics
            stats = cleanupTrainingStats(this);
            
            % create training result struct for analysis and recreating plot
            trainingInfoStruct = createTrainingInfoStruct(this,episodeIndex,trainingStoppedReason,trainingStoppedValue);
            stats.Information = trainingInfoStruct;
        end
        
        function stats = cleanupTrainingStats(this)
            maxEpisodes = this.TrainingOptions.MaxEpisodes;
            episodeIndex = this.EpisodeCount;
            % Clean up unused training statistics
            rmidx = (episodeIndex+1):maxEpisodes;
            stats = this.TrainingStats;
            if ~isempty(stats)
                stats.EpisodeIndex   (rmidx) = [];
                stats.EpisodeReward  (rmidx) = [];
                stats.EpisodeSteps   (rmidx) = [];
                stats.AverageReward  (rmidx) = [];
                stats.TotalAgentSteps(rmidx) = [];
                stats.AverageSteps   (rmidx) = [];
                if this.HasCritic
                    stats.EpisodeQ0  (rmidx) = [];
                end
                
                % attach the simulation info to the output structure
                stats.SimulationInfo = this.SimulationInfo;
                stats.SimulationInfo((episodeIndex+1):end) = [];
            end
            stats.SimulationInfo = vertcat(this.SimulationInfo{:});
        end
        function trainingInfoStruct = createTrainingInfoStruct(this,episodeIndex,trainingStoppedReason,trainingStoppedValue)
            % create training result strcut for analysis and plots
            if ~isempty(this.EpisodeMgr) && isvalid(this.EpisodeMgr) && this.EpisodeMgr.StopTrainingFromButtonFlag
                trainingStoppedReason = getString(message('rl:general:TextStopButton'));
                trainingStoppedValue = episodeIndex;
            end
            trainingInfoStruct = struct(...
                'TrainingOpts',                 this.TrainingOptions,...
                'HasCritic',                    this.HasCritic,...
                'HardwareResource',             this.Devices,...
                'LearningRate',                 this.LearnRates,...
                'TraningStartTime',             this.TrainingStartTime,...
                'ElapsedTime',                  this.TrainingElapsedTime,...
                'StopTrainingCriteria',         trainingStoppedReason,...
                'StopTrainingValue',            trainingStoppedValue...
                );
        end
        
        function train(this)
            % train the agent
            
            % create the trainer
            trainer = rl.train.createTrainerFactory(this.Environment,this.Agent,this.TrainingOptions);
            
            % attach the trainer to the training manager
            attachTrainer(this,trainer);
            % on cleanup, detatch the trainer
            cln = onCleanup(@() detatchTrainer(this,trainer));
            % run the trainer
            run(trainer);
        end
        function attachTrainer(this,trainer)
            % attach the training manager to a trainer

            this.TasksRunningOnWorkersListener = addlistener(trainer,'TasksRunningOnWorkers',...
                @(src,ed) setActionMessage(this,getString(message('rl:general:TrainingManagerRunningTasksOnWorkers'))));
            this.TasksCleanedUpOnWorkersListener = addlistener(trainer,'TasksRunningOnWorkers',...
                @(src,ed) setActionMessage(this,getString(message('rl:general:TrainingManagerCleaningUpWorkers'))));
            this.TasksActionMessageReceived = addlistener(trainer,'ActionMessageReceived',...
                @(src,ed) setActionMessage(this,ed.Data));
            
            % set the update fcn here (listeners will drop events if not
            % marked as recursive)
            trainer.FinishedEpisodeFcn   = @(info) update(this,info);
        end
        function detatchTrainer(this,trainer)
            % detatch the trainer from the training manager
            delete(trainer);
            
            delete(this.TasksRunningOnWorkersListener);
            delete(this.TasksCleanedUpOnWorkersListener);
            delete(this.TasksActionMessageReceived);
        end
    end
    methods (Access = private)
        function attachEnvEpisodeFinishedListener(this)
            this.EpisodeFinishedListener = addlistener(this.Environment,'EpisodeFinished',...
                    @(src,ed) update(this,ed.Data));
        end
        function initializeEpisodeManager(this)
            try
                %% get device and learning rate
                % REVISIT: agent abstraction does NOT define
                % getAction/getCritic
                % actor
                agent = this.Agent;
                if ~isempty(agent.getActor)
                    actorOptions = agent.getActor.Options;
                    actorDevice = actorOptions.UseDevice;
                    actorLearnRate = actorOptions.LearnRate;
                end
                % critic
                if ~isempty(agent.getCritic)
                    criticOptions = agent.getCritic.Options;
                    criticDevice = criticOptions.UseDevice;
                    criticLearnRate = criticOptions.LearnRate;
                end
                % three cass: actor only, critic only, both actor and critic
                if ~isempty(agent.getActor) && ~isempty(agent.getCritic)
                    devices = struct('actorDevice',actorDevice,'criticDevice',criticDevice);
                    learnRates = [actorLearnRate,criticLearnRate];
                elseif ~isempty(agent.getActor) && isempty(agent.getCritic)
                    devices = struct('actorDevice',actorDevice);
                    learnRates = actorLearnRate;
                elseif ~isempty(agent.getCritic) && isempty(agent.getActor)
                    devices = struct('criticDevice',criticDevice);
                    learnRates = criticLearnRate;
                end
            catch
                learnRates = 1;
                devices = struct('criticDevice','unknown');
            end
            this.Devices = devices;
            this.LearnRates = learnRates;
            % build the episode manager
            if strcmp(this.TrainingOptions.Plots,'training-progress')
                envname = getNameForEpisodeManager(this.Environment);

                %% create the episode manager
                delete(this.StopTrainingListener);
                episodeMgr = rl.util.EpisodeManager(this.TrainingOptions,this.Devices,this.LearnRates,this.HasCritic);
                agentname = regexprep(class(this.Agent),'\w*\.','');
                if isempty(envname)
                    emtitle = getString(message('rl:general:TrainingManagerEpisodeManagerTitle1'));
                else
                    emtitle = getString(message('rl:general:TrainingManagerEpisodeManagerTitle2',envname,agentname));
                end
                setTextOnAxis(episodeMgr,emtitle,...
                    getString(message('rl:general:TrainingManagerEpisodeManagerXLabel')),...
                    getString(message('rl:general:TrainingManagerEpisodeManagerYLabel')));

                this.StopTrainingListener = addlistener(episodeMgr,...
                    'RequestToStopTraining','PostSet',...
                    @(src,ed) request2ManuallyTerminateCB(this,src,ed));

                % store the episode manager
                this.EpisodeMgr = episodeMgr;
                setStartTime(this.EpisodeMgr,this.TrainingStartTime);
            end
        end
        function info = computeMetrics(this,epinfo)
            % returns relevant training progress metrics as a struct info.
            %
            %  Info.AverageSteps   : Running average of number of steps per episode
            %  Info.AverageReward  : Running average of reward per episode
            %  Info.EpisodeReward  : Reward for current episode
            %  Info.GlobalStepCount: Total times the agent was invoked
            %  Info.EpisodeCount   : Total number of episodes the agent has trained for

            episodeIndex = this.EpisodeCount;
            episodeSteps = epinfo.StepsTaken;
            episodeReward = epinfo.CumulativeReward;
            totalStepCount = this.TotalEpisodeStepCount;
            q0 = epinfo.Q0;
            
            % circular buffer index for averaging window
            numScoresToAverage = this.TrainingOptions.ScoreAveragingWindowLength;
            idx = mod(episodeIndex-1,numScoresToAverage)+1;
            this.RewardsInAveragingWindow(idx) = episodeReward;
            this.StepsInAveragingWindow(idx) = episodeSteps;
            
            numScores = min(episodeIndex,numScoresToAverage);
            avgReward = sum(this.RewardsInAveragingWindow)/numScores;
            avgSteps = sum(this.StepsInAveragingWindow)/numScores;
            
            info.AverageSteps    = avgSteps;
            info.AverageReward   = avgReward;
            info.EpisodeReward   = episodeReward;
            info.GlobalStepCount = totalStepCount;
            info.EpisodeCount    = episodeIndex;
            info.EpisodeSteps    = episodeSteps;
            if this.HasCritic
                info.EpisodeQ0 = q0;
            end
        end
        function updateCommandLineDisply(this,info)
            % update the command line display
            
            if this.TrainingOptions.Verbose
                episodeIndex = info.EpisodeCount;
                
                MaxEpisodes     = this.TrainingOptions.MaxEpisodes;
                stepCount       = info.EpisodeSteps;
                globalStepCount = info.GlobalStepCount;
                
                str = sprintf('Episode: %3d/%3d | Episode Reward : %4.2f | Episode Steps: %4d | Avg Reward : %.2f | Step Count : %d', ...
                    episodeIndex,MaxEpisodes,info.EpisodeReward,stepCount,info.AverageReward,globalStepCount);
                if this.HasCritic
                    str = sprintf('%s | Episode Q0 : %4.2f',str,info.EpisodeQ0);
                end
                fprintf('%s\n',str);
            end
        end
        function updateEpisodeManager(this,info)
            % push the training data onto the episode manager if
            % available
            episodeMgr = this.EpisodeMgr;
            if ~isempty(episodeMgr) && isvalid(episodeMgr) && ~this.Agent.TerminateSimulation
                stepEpisode(episodeMgr,info);               
            end
        end
        function updateTrainingStats(this,info)
            % Keep track of statistics
            stepCount       = info.EpisodeSteps;
            globalStepCount = info.GlobalStepCount;
            episodeIndex    = info.EpisodeCount;
            
            this.TrainingStats.EpisodeIndex   (episodeIndex) = episodeIndex;
            this.TrainingStats.EpisodeReward  (episodeIndex) = info.EpisodeReward;
            this.TrainingStats.EpisodeSteps   (episodeIndex) = stepCount;
            this.TrainingStats.AverageReward  (episodeIndex) = info.AverageReward;
            this.TrainingStats.TotalAgentSteps(episodeIndex) = globalStepCount;
            this.TrainingStats.AverageSteps   (episodeIndex) = info.AverageSteps;
            if this.HasCritic
                this.TrainingStats.EpisodeQ0  (episodeIndex) = info.EpisodeQ0;
            end
        end
        function saveAgentToDisk(this,info)
            % save the agent to disk if the provided criteria has been met
            episodeIndex = info.EpisodeCount;
            if this.TrainingOptions.SaveAgentFunction(info)
                SavedAgentFileName = fullfile(this.TrainingOptions.SaveAgentDirectory,['Agent' num2str(episodeIndex) '.mat']);
                saved_agent = this.Agent;
                stats = cleanupTrainingStats(this);
                savedAgentResultStruct = createSavedAgentResultStruct(this,stats); %#ok<NASGU>
                % make sure the save agent is in sim mode
                wasMode = getStepMode(saved_agent);
                setStepMode(saved_agent,"sim");
                try
                    save(SavedAgentFileName,'saved_agent', 'savedAgentResultStruct'); %#ok<USENS>
                catch
                    % g1928023: We do not want to interrupt the training
                    % due to saving errors. Therefore a warning is thrown.
                    warning(message('rl:general:TrainingManagerUnableToSaveAgent',this.TrainingOptions.SaveAgentDirectory))
                end
                % change the mode back
                setStepMode(saved_agent,wasMode);
            end
        end
        
        function savedAgentResultStruct = createSavedAgentResultStruct(this,stats)
            % elapsed time
            elapsedTime =  char(getDurationSinceReset(this.Watch));
            % saved agent info
            savedAgentInfoStruct = struct(...
                'TrainingOpts',                 this.TrainingOptions,...
                'HasCritic',                    this.HasCritic,...
                'HardwareResource',             this.Devices,...
                'LearningRate',                 this.LearnRates,...
                'TraningStartTime',             this.TrainingStartTime,...
                'ElapsedTime',                  elapsedTime);
            % saved agent result
            savedAgentResultStruct = struct(...
                'TrainingStats',                stats,...
                'Information',                  savedAgentInfoStruct);
        end
        
        function stopTraining = checkStopTraining(this,info)
            % stop training if the provided criteria has been met
            stopTraining = false;
            % Stop training (by stopping criteria or manually requested)
            if this.TrainingOptions.StopTrainingFunction(info) || ...
                    this.Agent.TerminateSimulation || ...
                    this.EpisodeCount >= this.TrainingOptions.MaxEpisodes
                stopTraining = true;
            end
        end
        function stopTraining = updateDisplaysFromTrainingInfo(this,info)
            % update the user visible components
            updateCommandLineDisply(this,info);
            updateEpisodeManager(this,info);
            % update training stats
            updateTrainingStats(this,info);
            % save agent to disk if requested
            saveAgentToDisk(this,info);
            % stop training
            stopTraining = checkStopTraining(this,info);
        end
        function request2ManuallyTerminateCB(this,~,ed) 
            % callback to manually terminate training
            terminateSimulation(this.Agent,ed.AffectedObject.RequestToStopTraining);
            % tell the world that training has been manually stopped
            notify(this,'TrainingManuallyStopped');
            % flush the event queue to make sure listeners can respond
            % immediately
            drawnow();
        end
    end
    methods(Hidden)
        function mgr = getEpisodeManager(this)
            mgr = this.EpisodeMgr;
        end
        function qeSaveAgentToDisk(this,info)
            saveAgentToDisk(this,info);                             
        end
    end
end

%% local utility functions
function str = iDateTimeAsStringUsingDefaultLocalFormat(dt)
defaultFormat = datetime().Format;
dt.Format = defaultFormat;
str = char(dt);
end
