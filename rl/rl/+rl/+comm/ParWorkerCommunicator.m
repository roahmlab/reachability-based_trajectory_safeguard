classdef ParWorkerCommunicator < rl.comm.AbstractWorkerCommunicator
% PARWORKERCOMMUNICATOR

% Revised: 10-16-2018
% Copyright 2018 The MathWorks, Inc.
    properties (Access = private)
        WorkerNode_
    end
    methods 
        function this = ParWorkerCommunicator(workerNode)
            this.WorkerNode_ = workerNode;
        end
        function delete(this)
            delete(this.WorkerNode_);
        end
        function id = getID(this)
            id = getID(this.WorkerNode_);
        end
        function params = handshakeExperience(this,exp)
            params = handshake(this.WorkerNode_,"ParametersUpdated","ExperienceReceived",exp);
        end
        function params = handshakeGradient(this,gradientBuffer)
            params = handshake(this.WorkerNode_,"ParametersUpdated","GradientReceived",gradientBuffer);
        end
        function params = handshakeParameterRequest(this)
            params = handshake(this.WorkerNode_,"ParametersUpdated","RequestedParameters",[]);
        end
        function action = handshakeRemoteGetAction(this,obs)
            action = handshake(this.WorkerNode_,"GetActionComputed","RemoteObservationsReceived",obs);
        end
        function action = handshakeRemoteGetInitialAction(this,obs)
            action = handshake(this.WorkerNode_,"GetInitialActionComputed","RemoteInitialObservationsReceived",obs);
        end
        function action = handshakeRemoteStep(this,exp)
            action = handshake(this.WorkerNode_,"StepActionComputed","RemoteExperienceReceived",exp);
        end
        function params = receiveParameters(this)
            % get the latest params from the communication channel without
            % blocking. Return [] if not in the buffer.
            
            % get the latest params from the communication channel without
            % blocking
            flag = true;
            params = [];
            while flag
                % get parameters back from the host
                p_ = getTopicData(this.WorkerNode_,"ParametersUpdated",false);
                if isempty(p_)
                    flag = false;
                else
                    params = p_;
                end
            end
        end
        function continueTraining = handshakeTrainingInfo(this,info)
            continueTraining = handshake(this.WorkerNode_,"ContinueTraining","SimulationFinished",info);
        end
        function sendTrainingInfo(this,info)
            sendData2Host(this.WorkerNode_,"SimulationFinished",info);
        end
        function cont = receiveContinueTraining(this,varargin)
            cont = getTopicData(this.WorkerNode_,"ContinueTraining",varargin{:});
        end
        function sendExperience(this,exp)
            sendData2Host(this.WorkerNode_,"ExperienceReceived",exp);
        end
        function sendGradient(this,gradientBuffer)
            sendData2Host(this.WorkerNode_,"GradientReceived",gradientBuffer);
        end
        function registerParametersReceivedCallback(this,agent,cb)
            addTopicCallback(this.WorkerNode_,"ParametersUpdated",@(msg) cb(this,agent,msg));
        end
        function registerContinueTrainingCallback(this,agent,cb)
            addTopicCallback(this.WorkerNode_,"ContinueTraining",@(msg) cb(this,agent,msg));
        end
        function registerManuallyStoppedTrainingCallback(this,agent,cb)
            addTopicCallback(this.WorkerNode_,"ManuallyStoppedTraining",@(msg) cb(this,agent,msg));
        end
    end
end