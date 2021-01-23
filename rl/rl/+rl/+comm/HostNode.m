classdef HostNode < rl.comm.CommNode
% HOSTNODE

% Revised: 10-9-2018
% Copyright 2017-2018 The MathWorks Inc.
    properties (Access = private)
        HostQueue
        WorkerQueues = []
        WorkerIDs = []
    end
    events
        WorkerRegistered
    end
    methods
        function this = HostNode()
            % build the host queue
            this = this@rl.comm.CommNode();
            this.HostQueue = parallel.pool.DataQueue();
            setupReceiveCB(this);
            
            % add a callback to populate the worker queues and IDs
            addTopicCallback(this,"InitializeWorker",@(msg)initializeWorkerTopicCB(this,msg));
        end
        function ids = getWorkerIDs(this)
            ids = this.WorkerIDs;
        end
        function hasAllIDs = blockUntilAllWorkersAreRegistered(this,workerIDs)
            % block until we receive messages from each of the workers
            % specified by workerIDs
            hasAllIDs = false;
            % add a timer exit flag to kill the blocking loop after X
            % seconds
            timerExit = false;
            t = timer();
            function nestedTimerFcn()
                timerExit = true;
            end
            t.TimerFcn = @(~,~) nestedTimerFcn();
            % REVISIT what is a good timeout to use here?
            TIMEOUT = 100;
            t.StartDelay = TIMEOUT;
            start(t);
            while (~hasAllIDs && ~timerExit)
                drawnow();
                hasAllIDs = areAllWorkersRegistered(this,workerIDs);
            end
            stop(t);
            delete(t);
        end
        function hasAllIDs = areAllWorkersRegistered(this,workerIDs)
            hasAllIDs = isempty(setdiff(workerIDs,getWorkerIDs(this)));
        end
        function hasAllIDs = areNWorkersRegistered(this,N)
            hasAllIDs = N == numel(getWorkerIDs(this));
        end
        function workerNode = createWorkerNode(this)
            % create a worker node, should be created by the execution of
            % the worker (in parfeval)
            workerNode = rl.comm.WorkerNode(this.HostQueue);
        end
        function sendData2Worker(this,topic,data,workerID)
            % send data to a worker specified by worker ID
            idx = getWorkerIdx(this,workerID);
            q = this.WorkerQueues(idx);
            msg = createMessage(this,topic,data);
            send(q,msg);
        end
        function sendData2AllWorkers(this,topic,data)
            % send data to all workers
            msg = createMessage(this,topic,data);
            for q = this.WorkerQueues(:)'
                send(q,msg);
            end
        end
    end
    methods (Access = protected)
        function queue = getIncomingQueue(this)
            queue = this.HostQueue;
        end
        function id = getWorkerID(~)
            % id = 0 -> maps to host, id > 0 maps to worker
            id = rl.util.getWorkerID();
        end
    end
    methods (Access = private)
        function idx = getWorkerIdx(this,id)
            idx = ismember(this.WorkerIDs,id);
        end
        function initializeWorkerTopicCB(this,msg)
            q  = msg.Data;
            id = msg.WorkerID;
            
            if any(ismember(this.WorkerIDs,id))
                error(message('rl:general:CommHostNodeDuplicateID',id));
            end
            
            this.WorkerQueues = [this.WorkerQueues;q ];
            this.WorkerIDs    = [this.WorkerIDs   ;id];
            
            % fire an event when a worker is registered
            ed = rl.util.RLEventData(id);
            notify(this,'WorkerRegistered',ed);
        end
    end
end