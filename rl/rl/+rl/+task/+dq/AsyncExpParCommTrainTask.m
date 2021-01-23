classdef AsyncExpParCommTrainTask < rl.task.dq.ParCommTrainTask
% AsyncExpParCommTrainTask
%
% This task sends async experiences to the host and expects parameters back
% via data queue communication objects

% Revised: 7-2-2019
% Copyright 2019 The MathWorks, Inc.

    properties (Access = private)
        % logger to build experience buffers
        ExpLogger
    end
    methods
        function this = AsyncExpParCommTrainTask(env,agent,trainOpts,workerComm)
            this = this@rl.task.dq.ParCommTrainTask(env,agent,trainOpts,workerComm);
            
            this.ExpLogger = rl.util.ExperienceLogger(this.StepsUntilDataIsSent);
            
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
            end
        end
        function processExperience(this,agent,exp)
            
            if ~this.UseBlocking
                % small pause to process parameter rcv callback
                pause(0.001);
            end
            
            addExperience2Buffer(this.ExpLogger,{exp});
            if isFull(this.ExpLogger) || (this.SendDataAtEndOfEpisode && exp{5} > 0)
                
                experienceBuffer = getExperiences(this.ExpLogger);
                experienceBuffer = preprocessExperience(agent,experienceBuffer);
                
                if this.UseBlocking
                    % send the experience back to the host and get parameters back
                    params = handshakeExperience(this.WorkerComm,experienceBuffer);
                    
                    % install the parameters from the host
                    if ~isempty(params)
                        setLearnableParameters(agent,params);
                    end
                    
                else
                    % send experiences back to host. parametes will be installed
                    % via the registered paramsReceivedCB callback
                    sendExperience(this.WorkerComm,experienceBuffer);
                end
                
                % reset the logger
                reset(this.ExpLogger);
            end
        end
    end
end