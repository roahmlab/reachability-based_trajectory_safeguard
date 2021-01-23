classdef rlPGAgent < rl.agent.AbstractOnPolicyPGAgent
    % rlPGAgent: Implements policy gradient (REINFORCE) agent
    
    % Copyright 2017-2019 The MathWorks Inc.
    
    properties (Access = private)
        % Flag used for validation since PG option validation differs
        % between pre-post construction
        IsConstructed = false
    end
    
    properties (Constant, Access = private)
        % Version indicator for backward compatibility
        Version = 2
        
        % Default buffer Length to create exp buffer at construction
        DefaultBufferLength = 500
    end
    
    methods
        function this = rlPGAgent(Actor, Critic, Options)
            % Constructor
            this = this@rl.agent.AbstractOnPolicyPGAgent(Actor, Critic, Options);
            
            this.IsConstructed = true;
        end
        
        function this = setCritic(this, Critic)
            % Overload superclass to change agent option UseBaseline opt
            this = setCritic@rl.agent.AbstractOnPolicyPGAgent(this,Critic);
            this.AgentOptions.UseBaseline = this.HasCritic;
        end
    end
    methods(Access = protected)
        function action = learn(this,exp)
            % learn from the current set of experiences where
            % exp = {state,action,reward,nextstate,isdone}
            % Return the noisy action.
            
            % store experiences
            appendExperience(this,exp);

            % accumulate gradient from experience
            gradientBuffer = accumulateGradient(this,[]);
            if ~isempty(gradientBuffer)
                % update parameters from learned gradients
                applyGradient(this,gradientBuffer)
                % reset experience buffer
                reset(this.ExperienceBuffer);
            end
            
            % compute action from the current policy
            % exp is {observation,action,reward,nextobservation,isdone}
            action = getActionWithExploration(this,exp{4});
        end
        
        function Actor = setActorLoss(this, Actor)
            if all(arrayfun(@(x) isa(x,'rl.util.rlNumericSpec'), this.ActionInfo))
                Actor = setLoss(Actor,"rl.loss.policyGradientContinuous");
            else
                Actor = setLoss(Actor,"rl.loss.policyGradientDiscrete");
            end
        end
        
        function validateOptionImpl(this,NewOptions)
            validateattributes(NewOptions,{'rl.option.rlPGAgentOptions'},{'scalar'},'','AgentOptions');
            if this.IsConstructed
                % NOTE: before construction, these checks happen on
                % informal API
                if NewOptions.UseBaseline && isempty(this.Critic)
                    error(message('rl:agent:errPGBaselineTrueNoCritic'))
                end
                if ~NewOptions.UseBaseline && ~isempty(this.Critic)
                    warning(message('rl:agent:warnPGBaselineFalseHaveCritic'))
                end
            end
        end
        
        function trainingOptions = validateAgentTrainingCompatibilityImpl(~,trainingOptions)
            % Validate PG agent parallel training options compatibility
            if ~strcmpi(trainingOptions.Parallelization,'none')
                dataToSend = trainingOptions.ParallelizationOptions.DataToSendFromWorkers;
                % PG agent only support send gradients for parallel
                if strcmpi(dataToSend,'Experiences')
                    error(message('rl:general:errParallelSendExpNotSupport'));
                end
                % validate StepsUntilDataIsSent
                if trainingOptions.ParallelizationOptions.StepsUntilDataIsSent ~= -1
                    error(message('rl:agent:errPGStepSendDataNotMinus1'));
                end
            end
        end
        
        function resetImpl(this)
            % Construct experience buffer. Buffer length is max train step
            % or simulation step
            
            if this.IsConstructed
                if (~isempty(this.Actor) && hasState(this.Actor)) || ...
                        (~isempty(this.Critic) && hasState(this.Critic))
                    error(message('rl:agent:errAgentHasStateNotSupport'))
                end
                
                % preferably build exp buffer with MaxSteps but when 
                % MaxStep is not built (resetImpl is called before train or
                % sim), build ExpBuffer with DefaultBufferLength
                if isinf(this.MaxSteps)
                    this.ExperienceBuffer = rl.util.ExperienceBuffer(this.DefaultBufferLength, this.ObservationInfo, this.ActionInfo);
                else
                    this.ExperienceBuffer = rl.util.ExperienceBuffer(this.MaxSteps, this.ObservationInfo, this.ActionInfo);
                end
                this.ExperienceBuffer.DoValidate = false;
            end
        end
    end
    
    methods (Hidden)
        function applyGradient(this,gradientBuffer)
            % Update representation from gradients
            
            if this.HasCritic
                this.Critic = optimize(this.Critic,gradientBuffer.CriticGrad);
            end
            this.Actor = optimize(this.Actor,gradientBuffer.ActorGrad);
        end
        
        function gradientBuffer = accumulateGradient(this,~)
            % Accumulate gradient from experience buffer
            % Return the a struct contains gradients of critic and actor.
            
            gamma = this.AgentOptions.DiscountFactor;
            % compute returns
            miniBatch = createExperienceWithReturnMiniBatch(this.ExperienceBuffer,gamma);
            if ~isempty(miniBatch)
                % compute gradients from sampled experiences
                gradientBuffer = stepRepresentation(this,miniBatch);
                % empty experience buffer at the end of the update interval
                reset(this.ExperienceBuffer);
            else
                gradientBuffer = [];
            end
        end
    end
    
    methods (Access = private)
        
        function gradientBuffer = stepRepresentation(this,MiniBatch)
            % Unpack experience
            Observations       = MiniBatch{1};
            Actions            = MiniBatch{2};
            ComputedReturns    = MiniBatch{3};
            TargetObservations = MiniBatch{4};
            IsDones            = MiniBatch{5};
            BatchSize  = numel(ComputedReturns);
            
            if this.AgentOptions.UseBaseline
                % Baseline target = n-step return + gamma^n*V[s+n]
                CurrentStateValues = getValue(this.Critic,Observations);
                TargetStateValue = getValue(this.Critic,TargetObservations);
                TargetStateValue(IsDones == 1) = 0; % early termination
                DiscountWeights = this.AgentOptions.DiscountFactor.^(BatchSize:-1:1);
                CriticTargets = ComputedReturns + DiscountWeights.*TargetStateValue;
                % Accumulate critic gradients in batch
                gradientBuffer.CriticGrad = gradient(this.Critic,'loss-parameters',Observations,CriticTargets);
                % Actor target = advantages = Baseline target - V[s]
                Advantages = reshape(CriticTargets,size(TargetStateValue)) - CurrentStateValues;
            else
                % Actor target = returns
                Advantages = ComputedReturns;
            end
            
            % Get actor targets for batch
            if isa(this.ActionInfo,'rl.util.rlNumericSpec')
                % continuous action
                % REVISIT: support multi continuous action channels
                LossVariable.Action = Actions{1};
                LossVariable.SamplingStrategy = this.Actor.SamplingStrategy;
            else
                % Actions indication matrix
                ActionIndicationMat = getElementIndicationMatrix(this.ActionInfo, Actions, BatchSize);
                % Actor targets
                Advantages = ActionIndicationMat .* Advantages;
            end
            
            % Accumulate actor gradients in batch
            LossVariable.Advantage = Advantages;
            LossVariable.EntropyLossWeight = this.AgentOptions.EntropyLossWeight;
            
            % Accumulate actor gradients in batch
            gradientBuffer.ActorGrad = gradient(this.Actor,'loss-parameters',Observations,LossVariable);
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
                    obj = rl.agent.rlPGAgent(Actor,Critic,s.AgentOptions_);
                end
            else
                obj = s;
            end
        end
    end
end
