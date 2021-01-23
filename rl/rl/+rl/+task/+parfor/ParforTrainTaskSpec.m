classdef ParforTrainTaskSpec < rl.task.TaskSpec
% PARFORTRAINTASKSPEC
%
% Specifications for parfor tasks

% Revised: 8-14-2019
% Copyright 2019 The MathWorks, Inc.
    properties
        TrainOpts
        Agent
        Env
    end
    properties (Access = private)
        Listeners = event.listener.empty
    end
    events
        WorkersSetup
        WorkersCleanedup
        EpisodeInfoReceivedFromWorker
        ActionMessageReceived
    end
    methods
        function this = ParforTrainTaskSpec(env,agent,trainOpts)
            this.Agent = agent;
            this.Env = env;
            this.TrainOpts = trainOpts;
        end
        function delete(this)
            delete@rl.task.TaskSpec(this);
            delete(this.Listeners)
        end
    end
    methods (Access = protected)
        % build the task
        function task = buildTask(this)
            popts = this.TrainOpts.ParallelizationOptions;
            if strcmpi(popts.DataToSendFromWorkers,'Gradients')
                task = rl.task.parfor.GradParforTrainTask(this.Env,this.Agent,this.TrainOpts);
            else
                task = rl.task.parfor.ExpParforTrainTask(this.Env,this.Agent,this.TrainOpts);
            end
            this.Listeners(1) = addlistener(task,"WorkersSetup"                 ,@(src,ed) notify(this,"WorkersSetup"                    ));
            this.Listeners(2) = addlistener(task,"WorkersCleanedup"             ,@(src,ed) notify(this,"WorkersCleanedup"                ));
            this.Listeners(3) = addlistener(task,"ActionMessageReceived"        ,@(src,ed) notify(this,"ActionMessageReceived"        ,ed));
            this.Listeners(4) = addlistener(task,"EpisodeInfoReceivedFromWorker",@(src,ed) notify(this,"EpisodeInfoReceivedFromWorker",ed));
        end
        % get the required number of workers for the task:
        %   0: no workers required (run on host)
        %   1: 1 worker required
        % Based on this number, the tasks will be ordered appropriately
        function n = getNumRequiredWorkers(~)
            n = 0;
        end
        % get the number of outputs for the task
        function n = getNumOutputs(~)
            n = 0;
        end
    end
end