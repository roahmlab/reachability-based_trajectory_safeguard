classdef ParCommTrainer < rl.train.Trainer
% ParCommTrainer
%
% Base class for data queue type trainers

% Revised: 7-2-2019
% Copyright 2019 The MathWorks, Inc.

    properties (Access = protected)
        % keep training or not
        ContinueTraining = true;
        % how many sims have finished
        FinishedSimCount = 0
        % how many sims have started
        StartedSimCount
        % flag to determine if the workers need to get parameters back as
        % soon as they are sent
        UseBlocking
        % count for the data received needed for sync
        ReceivedDataFromWorkerCount = 0
        % number of workers
        NumWorkers
    end
    methods (Abstract,Access = protected)
        tasks = buildTaskSpecs(this,hostComm,isAsync)
    end
    methods
        function this = ParCommTrainer(env,agent,trainOpts)
            this = this@rl.train.Trainer(env,agent,trainOpts);
            this.UseBlocking = this.TrainOpts.ParallelizationOptions.BlockWorkersUntilParametersAreReceived || ...
                ~strcmpi(this.TrainOpts.ParallelizationOptions.Mode,"async");
        end
        function run(this)
            % run the training
            pool = gcp();
            this.NumWorkers = pool.NumWorkers;
            this.StartedSimCount = this.NumWorkers;
            
            % build the hostComm
            hostComm = rl.comm.ParHostCommunicator();
            clnupHostComm = onCleanup(@() delete(hostComm));
            
            trainlist(1) = addlistener(this.Agent,'TerminateSimulation','PostSet',...
                @(src,ed) terminateSimCB(this,src,ed,hostComm));
            clnuplist = onCleanup(@() delete(trainlist));
            
            infcb = @(communicator,agent,message) ...
                trainingInfoReceivedFromWorkerCB(this,communicator,agent,message);
            registerTrainingInfoReceivedCallback(hostComm,this.Agent,infcb);
            
            isAsync = strcmpi(this.TrainOpts.ParallelizationOptions.Mode,"async");
            
            %% build the tasks
            taskSpecs = buildTaskSpecs(this,hostComm,isAsync);
            clntasks = onCleanup(@() delete(taskSpecs));
            
            %% run the tasks
            run(taskSpecs);
            
            % wait until N workers are registers
            while ~areNWorkersRegistered(hostComm,this.NumWorkers)
                pause(0.1);
            end
            sendContinueTraining(hostComm,true);
            
            notifyTasksRunning(this);
            clnup = onCleanup(@() notifyTasksCleanedUp(this));
            
            % wait for the tasks to finish
            waitForAllTasks2Finish(taskSpecs);
        end
    end
    methods (Access = protected)
        function terminateSimCB(~,~,ed,hostComm)
            if ed.AffectedObject.TerminateSimulation
                sendManuallyStopTraining(hostComm);
            end
        end
        function trainingInfoReceivedFromWorkerCB(this,communicator,agent,msg)
            % update the sim count
            this.FinishedSimCount = this.FinishedSimCount + 1;
            
            % get the worker id from the message
            wid = msg.WorkerID;
            
            % get the sim info from the worker
            info = msg.Data;
            info.WorkerID = wid;
            info.WorkerEpisodeCount = info.EpisodeCount;
            info.EpisodeCount = this.FinishedSimCount;
            
            % did the policy request to stop?
            this.ContinueTraining = ~agent.TerminateSimulation;
            
            % update the training manager
            if this.ContinueTraining && ...
                    (this.FinishedSimCount <= this.TrainOpts.MaxEpisodes)
                stopTraining = notifyEpisodeFinishedAndCheckStopTrain(this,info);
                this.ContinueTraining = ~stopTraining;
            else
                this.ContinueTraining = false;
            end
            if this.ContinueTraining && (this.StartedSimCount < this.TrainOpts.MaxEpisodes)
                sendContinueTraining(communicator,true,wid);
                this.StartedSimCount = this.StartedSimCount + 1;
            else
                sendContinueTraining(communicator,false,wid);
            end
        end
        function asyncSendParameters(this,communicator,agent,senderWID) %#ok<INUSL>
            % get the parameters from the agent
            p = getLearnableParameters(agent);

            % just send to the worker who sent the incoming data (sending
            % to all workers decreases sample efficiency by a bit)
            sendParameters(communicator,p,senderWID);
        end
        function syncSendParameters(this,communicator,agent)
            % record that data was received from the worker
            this.ReceivedDataFromWorkerCount = this.ReceivedDataFromWorkerCount + 1;
            if this.ReceivedDataFromWorkerCount >= getNumActiveSims(this)
                % NOTE parameters must be sent back to the workers even if
                % training was cancelled since they will block on the
                % worker until received
                
                % do some work before sending parameters (e.g. processing a
                % gradient buffer)
                syncPostReceiveDataFromAllWorkers(this,agent);
                
                % get the parameters from the agent
                p = getLearnableParameters(agent);
                
                % send the params back to all the workers
                sendParameters(communicator,p);
                
                % reset state
                this.ReceivedDataFromWorkerCount = 0;
            end
        end
        function syncPostReceiveDataFromAllWorkers(this,agent) %#ok<INUSD>
            % do work before sending parameters to simulation workers.
            % overloadable
        end
        function remainingEpisodes = getRemainingEpisodes(this)
            remainingEpisodes = this.TrainOpts.MaxEpisodes - this.FinishedSimCount;
        end
        function activeSims = getNumActiveSims(this)
            remainingEpisodes = getRemainingEpisodes(this);
            activeSims = min(this.NumWorkers,remainingEpisodes);
        end
    end
end