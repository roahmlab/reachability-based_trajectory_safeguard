classdef GradientBuffer < handle
% GRADIENTBUFFER
% Maintain a buffer of gradients to post-process (e.g. average)

% Revised: 10-2-2019
% Copyright 2019 The MathWorks, Inc.

    properties (Dependent,SetAccess = private)
        NumGradients
    end
    properties (SetAccess = private)
        % buffer of gradients
        GBuffer = []
    end
    methods
        function n = get.NumGradients(this)
            n = numel(this.GBuffer);
        end
        function append(this,g)
            % append the gradients to the buffer
            this.GBuffer = [this.GBuffer;g(:)];
        end
        function flush(this)
            % flush the buffer
            this.GBuffer = [];
        end
        function gavg = average(this)
            % average the elements in the bufffer. Currently supports
            % gradients in structure format (g.Critic, g.Actor)
            gavg = rl.internal.dataTransformation.averageGradients(this.GBuffer);
        end
    end
end