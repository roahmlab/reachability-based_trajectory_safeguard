classdef CommNode < handle
% COMMNODE

% Revised: 10-9-2018
% Copyright 2017-2018 The MathWorks Inc.

    properties (Access = private,Transient)
        % NOTE these properties MUST be transmit to prevent data from being
        % copied onto workers
        TopicBuffer(1,1) rl.comm.TopicBuffer
        TopicCallbacks
    end
    methods (Access = protected,Abstract)
        queue = getIncomingQueue(this)
        id = getWorkerID(this)
    end
    methods
        function this = CommNode()
            this.TopicBuffer = rl.comm.TopicBuffer();
            this.TopicCallbacks = containers.Map('KeyType','char','ValueType','any');
        end
        function [data,workerID] = getTopicData(this,topic,varargin)
            % get the first data and workerID that match the specified
            % topic from the TopicBuffer
            [data,workerID] = getTopicData(this.TopicBuffer,topic,varargin{:});
        end
        function addTopicCallback(this,topic,fcn)
            % add a callback to be fired once a message comes in with a
            % matching topic
            if isCallbackRegistered(this,topic)
                warning(message('rl:general:CommNodeOverwriteTopic',topic));
            end
            this.TopicCallbacks(topic) = fcn;
        end
        function rmTopicCallback(this,topic)
            % add a callback to be fired once a message comes in with a
            % matching topic
            if isCallbackRegistered(this,topic)
                remove(this.TopicCallbacks,topic);
            end
        end
        function val = isCallbackRegistered(this,topic)
            val = isKey(this.TopicCallbacks,topic);
        end
        function id = getID(this)
            id = getWorkerID(this);
        end
        function L = getQueueLength(this)
            q = getIncomingQueue(this);
            L = q.QueueLength;
        end
    end
    methods (Access = protected)
        function msg = createMessage(this,topic,data)
            % helper for creating outgoing messages
            msg = rl.comm.TopicMessage(topic,data,getWorkerID(this));
        end
        function setupReceiveCB(this)
            % setup receive callbacks on the queue which will receive data
            queue = getIncomingQueue(this);
            afterEach(queue,@(msg) receiveCB(this,msg))
        end
        function receiveCB(this,msg)
            % check to see if any callbacks with matching topics have been
            % added
            topic = msg.Topic;
            if isKey(this.TopicCallbacks,topic)
                cb = this.TopicCallbacks(topic);
                % fire the cb
                cb(msg);
            else
                % add the message to the buffer
                addMessage2Buffer(this.TopicBuffer,msg);
            end
            % FOR DEBUGGING
            % fprintf('the queue id = %d has %u elements in the buffer\n',getWorkerID(this),getQueueLength(this));
            % flush the event queue to prevent the queue from growing too
            % fast
            drawnow();
        end
    end

end