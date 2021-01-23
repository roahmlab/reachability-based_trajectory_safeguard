classdef ParCommTrainTask < rl.task.Task
% ParCommTrainTask

% Revised: 7-2-2019
% Copyright 2019 The MathWorks, Inc.

    properties
        TrainOpts
        Agent
        Env
        
        WorkerComm
        
        UseBlocking
        StepsUntilDataIsSent
        SendDataAtEndOfEpisode
    end
    methods (Abstract,Access = protected)
        processExperience(this,agent,exp)
    end
    methods
        function this = ParCommTrainTask(env,agent,trainOpts,workerComm)
            this = this@rl.task.Task();
            
            this.Agent = agent;
            this.Env = env;
            this.TrainOpts = trainOpts;
            % workerComm isn't built here for testing purposes
            this.WorkerComm = workerComm;
            
            popts                = trainOpts.ParallelizationOptions;
            stepsUntilDataIsSent = popts.StepsUntilDataIsSent;
            
            this.UseBlocking = popts.BlockWorkersUntilParametersAreReceived;
            if stepsUntilDataIsSent == -1
                this.StepsUntilDataIsSent = trainOpts.MaxStepsPerEpisode;
                this.SendDataAtEndOfEpisode = true;
            else
                this.StepsUntilDataIsSent = stepsUntilDataIsSent;
                this.SendDataAtEndOfEpisode = false;
            end
            
            % set the agent task functions
            setStepFcn(agent, ...
                @(agent,exp) taskStepFcn(this,agent,exp));
            setInitialActionFcn(agent, ...
                @(agent,obs) taskInitialActionFcn(this,agent,obs));
            % set the pre/post sim function on the agent
            setPreSimFcn(agent, ...
                @(agent,simCount) taskPreSimFcn(this,agent,simCount));
            setPostSimFcn(agent, ...
                @(agent,simCount,simInfo) taskPostSimFcn(this,agent,simCount,simInfo));
            
            % register a training stopped cb
            registerManuallyStoppedTrainingCallback(workerComm,agent,...
                @(comm,agent,msg) terminateSimCB(this,comm,agent,msg));
            
            drawnow();
            
        end
        function delete(this)
            delete@rl.task.Task(this);
            delete(this.WorkerComm);
            setStepMode(this.Agent,"sim");
        end
        function action = getActionWrapper(this,agent,obs) %#ok<INUSL>
            % overloadable
            action = getActionWithExploration(agent,obs);
        end
        function action = taskStepFcn(this,agent,exp)
            % step fcn

            % process the experience
            processExperience(this,agent,exp);
            
            % get action
            action = getActionWrapper(this,agent,exp{4});
        end
        function action = taskInitialActionFcn(this,agent,obs)
            % initial action fcn
            action = getActionWrapper(this,agent,obs);
        end
        function taskPreSimFcn(this,agent,simCount)
            % don't start simulating until a connection is established with the host
            if simCount < 2
                continueTraining = receiveContinueTraining(this.WorkerComm);
                terminateSimulation(agent,~continueTraining);
            end
        end
        function taskPostSimFcn(this,agent,simCount,simInfo)
            % transmit at the end of sim
            info.EpisodeInfo    = getEpisodeInfo(agent);
            info.EpisodeCount   = simCount;
            info.SimulationInfo = simInfo;
            continueTraining = handshakeTrainingInfo(this.WorkerComm,info);
            terminateSimulation(agent,~continueTraining);
        end
    end
    methods (Access = protected)
        function terminateSimCB(this,comm,agent,msg) %#ok<INUSD,INUSL>
            % DEBUG
            % save(sprintf('term_sim_%u',this.WorkerID),'msg');
            terminateSimulation(agent,true);
        end
        function varargout = runImpl(this)
            % run the simulations
            simOpts = getSimulationOptions(this.TrainOpts);
            simOpts.UseParallel = false;
            [varargout{1},varargout{2}] = simWithPolicy(this.Env,this.Agent,simOpts);
        end
    end
end