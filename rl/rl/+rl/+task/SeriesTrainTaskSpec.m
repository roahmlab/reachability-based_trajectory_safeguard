classdef SeriesTrainTaskSpec < rl.task.TaskSpec
% SIMPUBSUBTASKSPEC
%
% Specification to train an agent in series

% Revised: 6-25-2019
% Copyright 2019 The MathWorks, Inc.
    properties
        TrainOpts
        Agent
        Env
    end
    methods
        function this = SeriesTrainTaskSpec(env,agent,trainOpts)
            this.Agent = agent;
            this.Env = env;
            this.TrainOpts = trainOpts;
        end
    end
    methods (Access = protected)
        % build the task
        function task = buildTask(this)
            task = rl.task.SeriesTrainTask(this.Env,this.Agent,this.TrainOpts);
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
            n = 2;
        end
    end
end