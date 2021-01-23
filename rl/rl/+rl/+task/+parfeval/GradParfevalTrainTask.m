classdef GradParfevalTrainTask < rl.task.parfeval.ParfevalTrainTask
% GRADPARFORTRAINTASK
%
% Train an agent using parfeval, where each simulation sends gradients
% for learning. Parameters are installed at the top of
% every loop.

% Revised: 8-14-2019
% Copyright 2019 The MathWorks, Inc.

    properties (Access = private)
        GradientBuffer = rl.util.GradientBuffer();
    end
    methods
        function this = GradParfevalTrainTask(env,agent,trainOpts)
            this = this@rl.task.parfeval.ParfevalTrainTask(env,agent,trainOpts);
        end
    end
    methods (Access = protected)
        function agent = buildDistributedAgent(this)
            agent = this.Agent;
        end
        function setupParfevalWorkers(this)
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
        function learnFromWorkerDataImpl(this,wdata,requestedActiveSims)
            if isempty(requestedActiveSims)
                % apply the gradient right away
                applyGradient(this.Agent,wdata);
            else
                % append the gradiend to the buffer
                append(this.GradientBuffer,wdata);
                % if the buffer has requestedActiveSims gradients in the
                % buffer, average them and apply the gradient to the agent
                if this.GradientBuffer.NumGradients >= numSims
                    gavg = average(this.GradientBuffer);
                    applyGradient(this.Agent,gavg);
                    % flush the buffer
                    flush(this.GradientBuffer)
                end
            end
        end
        function str = generateDataReceivedMsgImpl(~,wdata,wid)
            str = getString(message(...
                'rl:general:TrainingManagerReceivedExperiences',...
                numel(wdata{1}),wid));
        end
    end
end