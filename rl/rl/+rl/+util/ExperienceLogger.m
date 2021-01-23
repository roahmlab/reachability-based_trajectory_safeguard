classdef ExperienceLogger < matlab.mixin.Copyable
% EXPERIENCELOGGER

% Copyright 2018 The MathWorks, Inc.

    properties (Access = private)
        Buffer
        BufferSize
        CurrentIndex
    end
    methods
        function this = ExperienceLogger(bufferSize)
            this.BufferSize = bufferSize;
            reset(this);
        end
        function val = isFull(this)
            val = (this.CurrentIndex - 1) >= this.BufferSize;
        end
        function experiences = getExperiences(this)
            % return the experiences as a cell
            experiences = this.Buffer(1:(this.CurrentIndex-1));
        end
        function reset(this)
            this.CurrentIndex = 1;
            this.Buffer = cell(this.BufferSize,1);
        end
        function addExperience2Buffer(this,experiences)
            % add experiences to the buffer
            szexp = numel(experiences);
            i1 = this.CurrentIndex;
            i2 = i1 + szexp - 1;
            if i2 > this.BufferSize
                error(message('rl:general:errExperienceLoggerNumElements'));
            end
            this.CurrentIndex = i2 + 1;
            this.Buffer(i1:i2) = experiences;
        end
    end
    methods (Access = private)
    end
end