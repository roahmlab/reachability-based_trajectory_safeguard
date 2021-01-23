classdef GradParCommTrainer < rl.train.dq.ParCommTrainer
% GradParCommTrainer
%
% Train an agent by sending gradients to the host via data queues.
% Parameters are sent back to the workers from the host. 

% Revised: 7-2-2019
% Copyright 2019 The MathWorks, Inc.
    
    properties (Access = private)
        GradientBuffer = rl.util.GradientBuffer();
    end
    methods
        function this = GradParCommTrainer(env,agent,trainOpts)
            this = this@rl.train.dq.ParCommTrainer(env,agent,trainOpts);
        end
    end
    methods (Access = protected)
        function asyncGradientsReceivedFromWorkerCB(this,communicator,agent,msg)
            % get the worker id from the message
            wid = msg.WorkerID;
            
            if ~agent.TerminateSimulation
                
                % get the gradients from the message
                gradients = msg.Data;
                
                % update representations from gradients (gorila)
                applyGradient(agent,gradients);
                
                % tell the world data is received
                str = getString(message(...
                    'rl:general:TrainingManagerReceivedGradients',...
                    wid));
                notifyActionMessageReceived(this,str);
            end
            
            % send params back to the worker
            asyncSendParameters(this,communicator,agent,wid);
        end
        function syncGradientsReceivedFromWorkerCB(this,communicator,agent,msg)
            % get the worker id from the message
            wid = msg.WorkerID;
            
            if ~agent.TerminateSimulation
                
                % get the gradients from the message
                gradients = msg.Data;
                
                % append the gradient to the buffer
                append(this.GradientBuffer,gradients);
                
                % tell the world data is received
                str = getString(message(...
                    'rl:general:TrainingManagerReceivedGradients',...
                    wid));
                notifyActionMessageReceived(this,str);
            end
            
            % send the params
            syncSendParameters(this,communicator,agent);
        end
        function syncPostReceiveDataFromAllWorkers(this,agent)
            % average the gradients
            gavg = average(this.GradientBuffer);
            
            % update representations from average gradients
            applyGradient(agent,gavg);
            
            % flush the buffer
            flush(this.GradientBuffer);
        end
        function tasks = buildTaskSpecs(this,hostComm,isAsync)
            % build the async exp task specs
            
            % register the experience received cb
            if isAsync
                datacb = @(communicator,agent,message) ...
                    asyncGradientsReceivedFromWorkerCB(this,communicator,agent,message);
            else
                datacb = @(communicator,agent,message) ...
                    syncGradientsReceivedFromWorkerCB(this,communicator,agent,message);
            end
            registerGradientReceivedCallback(hostComm,this.Agent,datacb);
            
            % build the task specs
            for i = 1:(this.NumWorkers)
                tasks(i) = rl.task.dq.AsyncGradParCommTrainTaskSpec(this.Env,this.Agent,this.TrainOpts,hostComm); %#ok<AGROW>
            end
        end
    end
end