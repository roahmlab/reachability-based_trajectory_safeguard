classdef WorkerNode < rl.comm.CommNode
% WORKERNODE

% Revised: 10-9-2018
% Copyright 2017-2018 The MathWorks Inc.

    properties (Access = private)
        HostQueue
        WorkerQueue
        WorkerID
    end
    methods
        function this = WorkerNode(hostQueue)
            % build the worker queue and attach the host queue
            this = this@rl.comm.CommNode();
            this.HostQueue = hostQueue;
            this.WorkerQueue = parallel.pool.DataQueue();
            % use current task id for worker ID for now
            this.WorkerID = rl.util.getWorkerID();
            % assign -1 for worker id if the current task is not in
            % parallel for testing purposes
            if this.WorkerID == 0
                this.WorkerID = -1;
            end
            setupReceiveCB(this);
            
            % send the created queue back to the host on the
            % "InitializeWorker" topic
            sendData2Host(this,"InitializeWorker",this.WorkerQueue);
        end
        function sendData2Host(this,topic,data)
            % send data back to the host node
            msg = createMessage(this,topic,data);
            send(this.HostQueue,msg);
        end
        function [data2recv,workerID] = handshake(this,recvtopic,sendtopic,data2send)
            % convenience function for establishing a blocking "handshake"
            % worker -> host -> worker
            sendData2Host(this,sendtopic,data2send);
            [data2recv,workerID] = getTopicData(this,recvtopic,true);
        end
    end
    methods (Access = protected)
        function queue = getIncomingQueue(this)
            queue = this.WorkerQueue;
        end
        function id = getWorkerID(this)
            id = this.WorkerID;
        end
    end
end