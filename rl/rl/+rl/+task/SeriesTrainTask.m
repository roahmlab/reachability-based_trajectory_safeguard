classdef SeriesTrainTask < rl.task.Task
% SERIESTRAINTASK
%
% Task to train an agent in series

% Revised: 7-2-2019
% Copyright 2019 The MathWorks, Inc.

    properties
        TrainOpts
        Agent
        Env
    end
    methods
        function this = SeriesTrainTask(env,agent,trainOpts)
            this = this@rl.task.Task();
            
            this.Agent = agent;
            this.Env = env;
            this.TrainOpts = trainOpts;
            
            setStepMode(agent,"learn");
        end
        function delete(this)
            delete@rl.task.Task(this);
            setStepMode(this.Agent,"sim");
        end
    end
    methods (Access = protected)
        function varargout = runImpl(this)
            simOpts = getSimulationOptions(this.TrainOpts);
            % return experiences and simulation outputs
            [varargout{1},varargout{2}] = simWithPolicy(this.Env,this.Agent,simOpts);
        end
    end
end