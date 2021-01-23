classdef TopicMessage
% TOPICMESSAGE

% Revised: 10-9-2018
% Copyright 2017-2018 The MathWorks Inc.

    properties
        % generic data sent
        Data
        % topic/channel for the data
        Topic string
        % worker ID should be an integer representing the worker the
        % message is being sent from. The host will be represented by a
        % WorkerID = 0
        WorkerID 
    end
    methods
        function this = TopicMessage(topic,data,workerID)
            this.Topic = topic;
            this.Data = data;
            this.WorkerID = workerID;
        end
        function idx = findTopicIndices(this,topic)
            % find the index of an array of messages with matching topic
            if isempty(this)
                idx = [];
            else
                topics = [this.Topic];
                idx = find(ismember(topics,string(topic)));
            end
        end
    end
end