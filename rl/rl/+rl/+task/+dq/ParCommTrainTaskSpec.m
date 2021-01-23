classdef ParCommTrainTaskSpec < rl.task.TaskSpec
% ParCommTrainTaskSpec

% Revised: 6-25-2019
% Copyright 2019 The MathWorks, Inc.
    properties
        TrainOpts
        Agent
        Env
        HostCommunicator
    end
    methods
        function this = ParCommTrainTaskSpec(env,agent,trainOpts,hostComm)
            this.Agent = agent;
            this.Env = env;
            this.TrainOpts = trainOpts;
            this.HostCommunicator = hostComm;
            
            % setup env if first instance
            if rl.task.dq.ParCommTrainTaskSpec.updateInstanceCount(1) == 0
                % REVISIT: this will move the simulations to temporary
                % directories in the simulink case.
                setupForSim(env,getSimulationOptions(trainOpts));
            end
        end
        function delete(this)
            delete@rl.task.TaskSpec(this);
            % cleanup env if last instance
            if rl.task.dq.ParCommTrainTaskSpec.updateInstanceCount(-1) == 1
                cleanupForSim(this.Env);
            end
        end
    end
    methods (Access = protected)
        % get the required number of workers for the task:
        %   0: no workers required (run on host)
        %   1: 1 worker required
        % Based on this number, the tasks will be ordered appropriately
        function n = getNumRequiredWorkers(~)
            n = 1;
        end
        % get the number of outputs for the task
        function n = getNumOutputs(~)
            n = 2;
        end
    end
    methods (Access = private,Static)
        function val = updateInstanceCount(d)
            persistent ct
            if isempty(ct)
                ct = 0;
            end
            val = ct;
            if nargin
                ct = ct + d;
            end
        end
    end
end