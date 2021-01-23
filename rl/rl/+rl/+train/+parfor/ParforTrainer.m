classdef ParforTrainer < rl.train.Trainer
% PARFORTRAINER
%
% Base class for parfor type trainers

% Revised: 8-14-2019
% Copyright 2019 The MathWorks, Inc.

    properties (Access = private)
        Listeners = event.listener.empty
    end
    methods
        function this = ParforTrainer(env,agent,trainOpts)
            this = this@rl.train.Trainer(env,agent,trainOpts);
        end
        function delete(this)
            delete(this.Listeners);
        end
        function run(this)
            
            % construct the task spec
            taskSpec = rl.task.parfor.ParforTrainTaskSpec(this.Env,this.Agent,this.TrainOpts);
            
            % attach listeners to the task
            this.Listeners(1) = addlistener(taskSpec,"WorkersSetup"                 ,@(src,ed) notifyTasksRunning  (this));
            this.Listeners(2) = addlistener(taskSpec,"WorkersCleanedup"             ,@(src,ed) notifyTasksCleanedUp(this));
            this.Listeners(3) = addlistener(taskSpec,"ActionMessageReceived"        ,@(src,ed) notifyActionMessageReceived           (this,ed));
            this.Listeners(4) = addlistener(taskSpec,"EpisodeInfoReceivedFromWorker",@(src,ed) notifyEpisodeFinishedAndCheckStopTrain(this,ed));
            
            % run the task
            run(taskSpec);
            
            % wait for the task to finish
            waitForAllTasks2Finish(taskSpec);
        end
    end
    methods (Access = private)
    end
end