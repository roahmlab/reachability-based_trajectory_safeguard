classdef AsyncGradParCommTrainTaskSpec < rl.task.dq.ParCommTrainTaskSpec
% AsyncGradParCommTrainTaskSpec
%
% Builds the Async Grads Task

% Revised: 6-25-2019
% Copyright 2019 The MathWorks, Inc.
    properties
    end
    methods
        function this = AsyncGradParCommTrainTaskSpec(env,agent,trainOpts,hostComm)
            this = this@rl.task.dq.ParCommTrainTaskSpec(env,agent,trainOpts,hostComm);
        end
    end
    methods (Access = protected)
        % build the task
        function task = buildTask(this)
            workerComm = buildParWorkerCommunicator(this.HostCommunicator);
            task = rl.task.dq.AsyncGradParCommTrainTask(this.Env,this.Agent,this.TrainOpts,workerComm);
        end
    end
end