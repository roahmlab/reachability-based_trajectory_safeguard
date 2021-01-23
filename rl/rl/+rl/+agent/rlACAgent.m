classdef rlACAgent < rl.agent.AbstractOnPolicyPGAgent
    % rlACAgent: Implements actor critic agent
    
    % Copyright 2017-2019 The MathWorks Inc.
    
    properties (Access = private)
        % Util to update step tracker
        LastUpdateStep
    end
    
    properties (Constant, Access = private)
        % version indicator for backward compatibility
        Version = 2
    end
    
    methods
        function this = rlACAgent(Actor, Critic, Options)
            % Constructor
            
            this = this@rl.agent.AbstractOnPolicyPGAgent(Actor, Critic, Options);
            
            % construct experience buffer and initialize step tracker
            resetImpl(this)
        end
    end
    
    methods(Access = protected)
        function Action = learn(this, Exp)
            % Learn from the current set of experiences where
            % exp = {state,action,reward,nextstate,isdone}
            % Return action with exploration.
            
            % store experiences
            appendExperience(this, Exp);
            
            StepCount = this.EpisodeInfo.StepsTaken;
            
            % update every n steps or if episode terminates
            BufferLength = StepCount - this.LastUpdateStep + 1;
            DoUpdate = (Exp{5} || BufferLength >= this.AgentOptions.NumStepsToLookAhead);
            if DoUpdate
                % compute gradients from n steps experiences in buffer
                % then reset buffer
                GradientBuffer = accumulateGradient(this, BufferLength);
                % apply gradient to actor and critic
                applyGradient(this, GradientBuffer);
                this.LastUpdateStep = StepCount + 1;
            end
            
            if Exp{5}
                % reset step tracker
                this.LastUpdateStep = 1;
            end
            
            % compute action from the current policy
            % exp is {observation,action,reward,nextobservation,isdone}
            Action = getActionWithExploration(this, Exp{4});
        end
        
        function Actor = setActorLoss(this, Actor)
            if all(arrayfun(@(x) isa(x,'rl.util.rlNumericSpec'), this.ActionInfo))
                Actor = setLoss(Actor,"rl.loss.policyGradientContinuous");
            else
                Actor = setLoss(Actor,"rl.loss.policyGradientDiscrete");
            end
        end
        
        function validateOptionImpl(~, NewOptions)
            validateattributes(NewOptions,{'rl.option.rlACAgentOptions'},{'scalar'},'','AgentOptions');
        end
        
        function trainingOptions = validateAgentTrainingCompatibilityImpl(this,trainingOptions)
            % Validate AC agent training options compatibility
            
            if ~strcmpi(trainingOptions.Parallelization,'none')
                dataToSend = trainingOptions.ParallelizationOptions.DataToSendFromWorkers;
                % AC agent only support send gradients for parallel
                if strcmpi(dataToSend,'Experiences')
                    error(message('rl:general:errParallelSendExpNotSupport'));
                end
                % validate StepsUntilDataIsSent
                if trainingOptions.ParallelizationOptions.StepsUntilDataIsSent == -1
                    error(message('rl:agent:errACStepSendDataFullEpi'));
                end
                % validate NumStepsToLookAhead and StepsUntilDataIsSent
                if trainingOptions.ParallelizationOptions.StepsUntilDataIsSent ~= this.AgentOptions.NumStepsToLookAhead
                    warning(message('rl:agent:warnACNumStepLookAheadNotUseInParallel'));
                end
            end
        end
        
        function resetImpl(this)
            % Construct experience buffer, should only be called after
            % construction (TrainingOptions only available when train is
            % called)
            
            if (~isempty(this.Actor) && hasState(this.Actor)) || ...
                    (~isempty(this.Critic) && hasState(this.Critic))
                error(message('rl:agent:errAgentHasStateNotSupport'))
            end
            
            this.ExperienceBuffer = rl.util.ExperienceBuffer(this.AgentOptions.NumStepsToLookAhead, this.ObservationInfo, this.ActionInfo);
            this.ExperienceBuffer.DoValidate = false;
            this.LastUpdateStep = 1;
        end
    end
    
    methods (Hidden)
        function applyGradient(this,gradientBuffer)
            % Update representation from gradients
            
            this.Critic = optimize(this.Critic,gradientBuffer.CriticGrad);
            this.Actor  = optimize(this.Actor,gradientBuffer.ActorGrad);
        end
        
        function GradientBuffer = accumulateGradient(this, BufferTime)
            % Compute gradient from the experience buffer
            % Return the a struct contains gradients of critic and actor.
            
            % Compute finite horizon advantage
            [Advantages, CriticTargets, BatchExperience] = computeFiniteHorizonAdvantage(this.ExperienceBuffer, ...
                        this.Critic, this.AgentOptions.DiscountFactor);
            
            % Unpack experience
            Observation = BatchExperience{1};
            Action      = BatchExperience{2};
            if isa(this.ActionInfo,'rl.util.rlNumericSpec')
                % continuous action
                % REVISIT: support multi continuous action channels
                LossVariable.Action = Action{1};
                LossVariable.SamplingStrategy = this.Actor.SamplingStrategy;
            else
                % Actions indication matrix
                ActionIndicationMat = getElementIndicationMatrix(this.ActionInfo, Action, BufferTime);
                % Actor targets
                Advantages = ActionIndicationMat .* reshape(Advantages,1,[]);
            end
            
            % Accumulate critic gradients in batch
            GradientBuffer.CriticGrad = gradient(this.Critic,'loss-parameters',Observation,CriticTargets);
            % Accumulate actor gradients in batch
            LossVariable.Advantage = Advantages;
            LossVariable.EntropyLossWeight = this.AgentOptions.EntropyLossWeight;
            GradientBuffer.ActorGrad = gradient(this.Actor,'loss-parameters',Observation,LossVariable);
            % empty experience buffer at the end of the update interval
            reset(this.ExperienceBuffer);
        end
    end
    
    methods (Static)
        function obj = loadobj(s)
            if isstruct(s)
                if ~isfield(s,'Version')
                    % version 1 does not have Version field
                    % In version 2,
                    %   - Critic changes from rlRepresentation to rlValueRepresentation
                    %   - Actor changes from rlRepresentation to rlStochasticActorRepresentation
                    %   - ExperienceBuffer requires obs and act dims inputs
                    %   but will always reconstruct since agent is on-policy
                    [Actor, Critic] = rl.agent.AbstractOnPolicyPGAgent.redirectV1Rep(s);
                    obj = rl.agent.rlACAgent(Actor,Critic,s.AgentOptions_);
                end
            else
                obj = s;
            end
        end
    end
end
