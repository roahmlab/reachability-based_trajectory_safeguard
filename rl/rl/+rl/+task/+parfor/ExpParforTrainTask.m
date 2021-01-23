classdef ExpParforTrainTask < rl.task.parfor.ParforTrainTask
% EXPPARFORTRAINTASK
%
% Train an agent using parfor, where each simulation send experiences back
% to the host for learning. Parameters are installed at the top of every
% parfor loop.

% Revised: 8-14-2019
% Copyright 2019 The MathWorks, Inc.

    methods
        function this = ExpParforTrainTask(env,agent,trainOpts)
            this = this@rl.task.parfor.ParforTrainTask(env,agent,trainOpts);
        end
    end
    methods (Access = protected)
        function agent = buildDistributedAgent(this)
            agent = copy(this.Agent);
            % REVISIT it would be good to put DDPG and DQN under a common
            % subclass to prevent having to check types
            if isa(agent.AgentOptions,'rl.option.AgentMemoryTarget')
                % reset the experience buffer as it has no impact on
                % training and can cause memory issues
                reset(agent.ExperienceBuffer);
            end
        end
        function setupParforWorkers(this)
            % run in exploration mode and attach a logger
            setStepMode(this.DistributedAgent,"sim-with-exploration");
            attachLogger(this.DistributedAgent,this.SimOpts.MaxSteps);
        end
        function [wdata,wsimInfo,wepInfo] = runWorkerSims(this)
            % sim
            [wdata,wsimInfo] = simWithPolicy(this.Env,this.DistributedAgent,this.SimOpts);
            wdata{1} = preprocessExperience(this.DistributedAgent,wdata{1});
            % get the episode info from the agent
            wepInfo = getEpisodeInfo(this.DistributedAgent);
        end
        function learnFromWorkerDataImpl(this,wdata,~)
            learnFromExperiences(this.Agent,wdata{1});
        end
        function str = generateDataReceivedMsgImpl(~,wdata,wid)
            str = getString(message(...
                'rl:general:TrainingManagerReceivedExperiences',...
                numel(wdata{1}),wid));
        end
    end
end