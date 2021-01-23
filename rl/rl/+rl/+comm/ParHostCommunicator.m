classdef ParHostCommunicator < rl.comm.AbstractHostCommunicator
% PARHOSTCOMMUNICATOR

% Revised: 10-16-2018
% Copyright 2018 The MathWorks, Inc.
    properties (Access = private)
        HostNode_ (1,1) rl.comm.HostNode
    end
    events
        WorkerRegistered
    end
    methods
        function this = ParHostCommunicator()
            this.HostNode_ = rl.comm.HostNode();
            addlistener(this.HostNode_,'WorkerRegistered',@(src,ed) workerRegisteredBridgeCB(this,src,ed));
        end
        function delete(this)
            delete(this.HostNode_);
        end
        function ids = getWorkerIDs(this)
            ids = getWorkerIDs(this.HostNode_);
        end
        function hasAllIDs = blockUntilAllWorkersAreRegistered(this,workerIDs)
            hasAllIDs = blockUntilAllWorkersAreRegistered(this.HostNode_,workerIDs);
        end
        function hasAllIDs = areAllWorkersRegistered(this,workerIDs)
            hasAllIDs = areAllWorkersRegistered(this.HostNode_,workerIDs);
        end
        function hasAllIDs = areNWorkersRegistered(this,N)
            hasAllIDs = areNWorkersRegistered(this.HostNode_,N);
        end
        function comm = buildParWorkerCommunicator(this)
            % SHOULD BE CREATED ON WORKER
            workerNode = createWorkerNode(this.HostNode_);
            comm = rl.comm.ParWorkerCommunicator(workerNode);
        end
        % callback syntax:
        %   callback(cb,comm,agent,message)
        function registerExperienceReceivedCallback(this,agent,cb)
            addTopicCallback(this.HostNode_,"ExperienceReceived",@(msg) cb(this,agent,msg));
        end
        function registerGradientReceivedCallback(this,agent,cb)
            addTopicCallback(this.HostNode_,"GradientReceived",@(msg) cb(this,agent,msg));
        end
        function registerTrainingInfoReceivedCallback(this,agent,cb)
            addTopicCallback(this.HostNode_,"SimulationFinished",@(msg) cb(this,agent,msg));
        end
        function registerParametersRequestedCallback(this,agent,cb)
            addTopicCallback(this.HostNode_,"RequestedParameters",@(msg) cb(this,agent,msg));
        end
        function registerRemoteGetActionCallback(this,agent,cb)
            addTopicCallback(this.HostNode_,"RemoteObservationsReceived",@(msg) cb(this,agent,msg));
        end
        function registerRemoteGetInitialActionCallback(this,agent,cb)
            addTopicCallback(this.HostNode_,"RemoteInitialObservationsReceived",@(msg) cb(this,agent,msg));
        end
        function registerRemoteStepCallback(this,agent,cb)
            addTopicCallback(this.HostNode_,"RemoteExperienceReceived",@(msg) cb(this,agent,msg));
        end
        function sendRemoteGetAction(this,action,workerID)
            sendData2Worker(this.HostNode_,"GetActionComputed",action,workerID);
        end
        function sendRemoteGetInitialAction(this,action,workerID)
            sendData2Worker(this.HostNode_,"GetInitialActionComputed",action,workerID);
        end
        function sendRemoteStep(this,action,workerID)
            sendData2Worker(this.HostNode_,"StepActionComputed",action,workerID);
        end
        function sendParameters(this,params,workerID)
            if nargin < 3
                sendData2AllWorkers(this.HostNode_,"ParametersUpdated",params);
            else
                sendData2Worker(this.HostNode_,"ParametersUpdated",params,workerID);
            end
        end
        function sendContinueTraining(this,val,workerID)
            if nargin < 3
                sendData2AllWorkers(this.HostNode_,"ContinueTraining",val);
            else
                sendData2Worker(this.HostNode_,"ContinueTraining",val,workerID);
            end
        end
        function sendManuallyStopTraining(this)
            sendData2AllWorkers(this.HostNode_,"ManuallyStoppedTraining",[]);
        end
    end
    methods (Access = private)
        function workerRegisteredBridgeCB(this,~,ed)
            notify(this,'WorkerRegistered',ed);
        end
    end
end