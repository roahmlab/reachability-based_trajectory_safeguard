classdef ExpParCommTrainer < rl.train.dq.ParCommTrainer
% ExpParCommTrainer
%
% Train an agent by sending experiences to the host via data queues.
% Parameters are sent back to the workers from the host. 

% Revised: 7-2-2019
% Copyright 2019 The MathWorks, Inc.

    methods
        function this = ExpParCommTrainer(env,agent,trainOpts)
            this = this@rl.train.dq.ParCommTrainer(env,agent,trainOpts);
        end
    end
    methods (Access = protected)
        function learnFromExperiences(this,agent,msg,wid)
            if ~agent.TerminateSimulation
                
                % get the experience from the message
                experiences = msg.Data;
                % learn from the experience
                learnFromExperiences(agent,experiences);
                % tell the world data is received
                str = getString(message(...
                    'rl:general:TrainingManagerReceivedExperiences',...
                    numel(experiences),wid));
                notifyActionMessageReceived(this,str);
            end
        end
        function asyncExperienceReceivedFromWorkerCB(this,communicator,agent,msg)
            % get the worker id from the message
            wid = msg.WorkerID;
            
            % learn from experiences
            learnFromExperiences(this,agent,msg,wid);
            
            % send params back to the worker
            asyncSendParameters(this,communicator,agent,wid);
        end
        function syncExperienceReceivedFromWorkerCB(this,communicator,agent,msg)
            % get the worker id from the message
            wid = msg.WorkerID;
            
            % learn from experiences
            learnFromExperiences(this,agent,msg,wid);
            
            % send the params
            syncSendParameters(this,communicator,agent);
        end
        function tasks = buildTaskSpecs(this,hostComm,isAsync)
            % build the async exp task specs
            
            % register the experience received cb
            if isAsync
                expcb = @(communicator,agent,msg) ...
                    asyncExperienceReceivedFromWorkerCB(this,communicator,agent,msg);
            else
                expcb = @(communicator,agent,msg) ...
                    syncExperienceReceivedFromWorkerCB(this,communicator,agent,msg);
            end
            registerExperienceReceivedCallback(hostComm,this.Agent,expcb);
            
            % make a copy of the agent to distribute
            cagent = copy(this.Agent);
            
            % REVISIT it would be good to put DDPG and DQN under a common
            % subclass to prevent having to check types
            if isa(cagent.AgentOptions,'rl.option.AgentMemoryTarget')
                % reset the experience buffer as it has no impact on
                % training and can cause memory issues
                reset(cagent.ExperienceBuffer);
            end
            
            % build the task specs
            for i = 1:(this.NumWorkers)
                tasks(i) = rl.task.dq.AsyncExpParCommTrainTaskSpec(this.Env,cagent,this.TrainOpts,hostComm); %#ok<AGROW>
            end
        end
    end
end