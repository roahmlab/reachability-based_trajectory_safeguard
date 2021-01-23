classdef GradParforTrainTask < rl.task.parfor.ParforTrainTask
% GRADPARFORTRAINTASK
%
% Train an agent using parfor, where each simulation send gradients back
% to the host for learning. Parameters are installed at the top of every
% parfor loop.

% Revised: 8-14-2019
% Copyright 2019 The MathWorks, Inc.

    properties
        GradientBuffer = rl.util.GradientBuffer()
    end
    methods
        function this = GradParforTrainTask(env,agent,trainOpts)
            this = this@rl.task.parfor.ParforTrainTask(env,agent,trainOpts);
        end
    end
    methods (Access = protected)
        function setupParforWorkers(this)
            % run in exploration mode and append experiences
            setStepMode(this.DistributedAgent,"sim-with-exploration-and-append");
        end
        function [wdata,wsimInfo,wepInfo] = runWorkerSims(this)
            % sim
            [~,wsimInfo] = simWithPolicy(this.Env,this.DistributedAgent,this.SimOpts);
            % get the episode info from the agent
            wepInfo = getEpisodeInfo(this.DistributedAgent);
            % accumulate the gradients from the learned episode
            wdata = accumulateGradient(this.DistributedAgent,wepInfo.StepsTaken);
        end
        function learnFromWorkerDataImpl(this,wdata,numSims)
            % append the gradiend to the buffer
            append(this.GradientBuffer,wdata);
            % if the buffer has numSims gradients in the buffer, average
            % them and apply the gradient to the agent
            if this.GradientBuffer.NumGradients >= numSims
                gavg = average(this.GradientBuffer);
                applyGradient(this.Agent,gavg);
                % flush the buffer
                flush(this.GradientBuffer)
            end
        end
        function str = generateDataReceivedMsgImpl(~,~,wid)
            str = getString(message(...
                'rl:general:TrainingManagerReceivedGradients',...
                wid));
        end
    end
end