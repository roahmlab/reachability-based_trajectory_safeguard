classdef AsyncGradParCommTrainTask < rl.task.dq.ParCommTrainTask
% AsyncGradParCommTrainTask
%
% This task sends async gradients to the host and expects parameters back
% via data queue communication objects

% Revised: 7-2-2019
% Copyright 2019 The MathWorks, Inc.

    properties (Access = private)
        % step count
        StepCount = 0
        % check to see if new params are received (initialized as true
        % since the latest params are available at the start of training)
        NewParametersReceived = true
    end
    methods
        function this = AsyncGradParCommTrainTask(env,agent,trainOpts,workerComm)
            this = this@rl.task.dq.ParCommTrainTask(env,agent,trainOpts,workerComm);
            
            if ~this.UseBlocking
                paramrcvcb = @(comm,agent,message) paramsReceivedCB(this,comm,agent,message);
                registerParametersReceivedCallback(this.WorkerComm,agent,paramrcvcb);
            end
        end
    end
    methods (Access = protected)
        function paramsReceivedCB(this,comm,agent,msg) %#ok<INUSL>
            % install the parameters sent by the host on the worker
            params = msg.Data;
            if ~isempty(params)
                setLearnableParameters(agent,params);
                this.NewParametersReceived = true;
            end
        end
        function processExperience(this,agent,exp)
            
            appendExperience(agent,exp);
            this.StepCount = this.StepCount + 1;
            if (this.StepCount >= this.StepsUntilDataIsSent) || (this.SendDataAtEndOfEpisode && exp{5} > 0)
                
                if ~this.UseBlocking
                    % wait until we have the latest params before gradients
                    % are computed
                    while ~this.NewParametersReceived
                        % small pause to process parameter rcv callback
                        pause(0.001);
                    end
                end

                % accumulate gradients from agent's experience buffer
                gradientBuffer = accumulateGradient(agent,this.StepCount);
                
                if ~isempty(gradientBuffer)
                    if this.UseBlocking
                        % send the gradient to the host and get parameters back
                        params = handshakeGradient(this.WorkerComm,gradientBuffer);

                        % install the parameters from the host
                        setLearnableParameters(agent,params);
                    else
                        % send gradients back to the host, parameters will come
                        % back via callback
                        sendGradient(this.WorkerComm,gradientBuffer);

                        % reset the new params rcv flag
                        this.NewParametersReceived = false;
                    end
                end
                
                % reset the step count and experience buffer
                this.StepCount = 0;
            end
        end
    end
end