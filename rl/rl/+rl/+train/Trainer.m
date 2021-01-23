classdef Trainer < handle
% TRAINER
%
% Base class for various training methods. Provides various utilities for
% communicating with the TrainingManager (see attachTrainer) via events and
% callbacks. Subclasses must define the run method.

% Revised: 7-2-2019
% Copyright 2019 The MathWorks, Inc.

    properties
        Env
        Agent
        TrainOpts
    end
    properties (Hidden)
        % updated in attachTrainer from the training manager
        FinishedEpisodeFcn = [] 
        
        % NOTE: the function handle approach vs the listeners seems like it
        % might be the right way to go as events seem to be dropped
        % especially in the pub sub case if the listeners aren't set to
        % recursive
    end
    events
        EpisodeFinished
        DataReceivedFromWorker
        TasksRunningOnWorkers
        TasksCleanedUpOnWorkers
        ActionMessageReceived
    end
    methods (Abstract)
        run(this)
    end
    methods 
        function this = Trainer(env,agent,trainOpts)
            this.Env = env;
            this.Agent = agent;
            this.TrainOpts = trainOpts;
        end
        function notifyEpisodeFinished(this,info)
            % fire an event when an episode is finished
            if isa(info,'rl.util.RLEventData')
                ed = info;
            else
                ed = rl.util.RLEventData(info);
            end
            notify(this,"EpisodeFinished",ed);
        end
        function stopTraining = notifyEpisodeFinishedAndCheckStopTrain(this,info)
            % fire an event when an episode is finished and check if
            % training should stop
            if isempty(this.FinishedEpisodeFcn)
                notifyEpisodeFinished(this,info);
                stopTraining = this.Agent.TerminateSimulation;
            else
                if isa(info,'rl.util.RLEventData')
                    info = info.Data;
                end
                stopTraining = this.FinishedEpisodeFcn(info);
            end
        end
        function notifyDataReceivedFromWorker(this,data)
            % fire an event when data is received from a worker. This
            % should only be used for testing (e.g. TestingMode = true)
            if isa(data,'rl.util.RLEventData')
                ed = data;
            else
                ed = rl.util.RLEventData(data);
            end
            notify(this,"DataReceivedFromWorker",ed);
        end
        function notifyTasksRunning(this)
            % fire an event telling the world tasks are running
            notify(this,"TasksRunningOnWorkers");
        end
        function notifyTasksCleanedUp(this)
            % fire an event telling the world tasks have been destroyed
            notify(this,"TasksCleanedUpOnWorkers");
        end
        function notifyActionMessageReceived(this,msg)
            % fire an event to redirect messages to the episode manager
            if isa(msg,'rl.util.RLEventData')
                ed = msg;
            else
                ed = rl.util.RLEventData(msg);
            end
            notify(this,"ActionMessageReceived",ed);
        end
    end

end