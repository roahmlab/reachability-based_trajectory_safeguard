classdef AsyncExpParCommTrainTaskSpec < rl.task.dq.ParCommTrainTaskSpec
% AsyncExpParCommTrainTaskSpec
%
% Builds the Async Exp Task

% Revised: 6-25-2019
% Copyright 2019 The MathWorks, Inc.
    properties
    end
    methods
        function this = AsyncExpParCommTrainTaskSpec(env,agent,trainOpts,hostComm)
            this = this@rl.task.dq.ParCommTrainTaskSpec(env,agent,trainOpts,hostComm);
        end
    end
    methods (Access = protected)
        % build the task
        function task = buildTask(this)
            workerComm = buildParWorkerCommunicator(this.HostCommunicator);
            task = rl.task.dq.AsyncExpParCommTrainTask(this.Env,this.Agent,this.TrainOpts,workerComm);
        end
    end
end