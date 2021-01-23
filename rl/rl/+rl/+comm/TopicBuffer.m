classdef TopicBuffer < handle
% TOPICBUFFER

% Revised: 10-9-2018
% Copyright 2017-2018 The MathWorks Inc.
    properties (Access = private)
        Buffer rl.comm.TopicMessage
    end
    properties (SetAccess = private,Dependent)
        BufferSize
    end
    events
        MessageReceived
    end
    methods
        function this = TopicBuffer()
            this.Buffer = rl.comm.TopicMessage.empty(1,0);
        end
        function val = get.BufferSize(this)
            val = numel(this.Buffer);
        end
        function b = getBuffer(this)
            b = this.Buffer;
        end
        function addMessage2Buffer(this,msg)
            this.Buffer = vertcat(this.Buffer,msg(:));
            notify(this,'MessageReceived');
        end
        function [data,workerID] = getTopicData(this,topic,isblocking)
            if nargin < 3
                isblocking = true;
            end
            
            idx = findTopicIndices(this.Buffer,topic);
            if isblocking
                while isempty(idx)
                    % REVISIT is drawnow slowing things down? In any case
                    % it shouldn't be necessary. 
                    % drawnow();
                    idx = findTopicIndices(this.Buffer,topic);
                end
            end
            if isempty(idx)
                data = [];
                workerID = [];
            else
                idx_ = idx(1);
                msg = this.Buffer(idx_);
                data = msg.Data;
                workerID = msg.WorkerID;
                
                % pop the message off the buffer
                this.Buffer(idx_) = [];
            end
        end
    end
end