classdef rlSACAgentNoTanh < rl.agent.AbstractAgent
    % rlSACAgentNoTanh: Implements a soft actor critic with
    % adaptive entropy weight
    % REVISIT: support RNN
    
    % Copyright 2019 The MathWorks Inc.
    
    properties (Dependent)
        % Options to configure RL agent
        AgentOptions
    end
    
    properties (Access = private)
        % Private options to configure RL agent
        AgentOptions_ = [];
        % Counter determined when to update the targets
        TargetUpdateCounter
    end
    
    properties (SetAccess = private)
        % Circular buffer for experience replay
        ExperienceBuffer
    end
    properties (Access = private)
        % Critic function approximators
        Critic1 
        Critic2
        
        % Target Critic function approximators
        TargetCritic1
        TargetCritic2
        
        % Actor function approximator
        Actor
        
        % Entropy related properties
        % Indicator determining tunable Entropy Weight
        LearnableEntropy
        
        % Current Entropy Weight
        EntropyWeight
       
        % optimize log of entropy weight to ensure non-negative values
        LogEntropyWeight
        
        TargetLogEntropy
        
        EntropyOptimizer
    end
    methods
        function this = rlSACAgentNoTanh(Actor, Critic1, varargin)
            % Constructor
            this = this@rl.agent.AbstractAgent();
            
            % extract observation and action info
            this.ActionInfo = Actor.ActionInfo;
            this.ObservationInfo = Actor.ObservationInfo;
            
            % REVISIT: Not support multiple action channels due to DLT
            % limitation (a loss layer cannot take multiple inputs)
            % multi actions might work with dlnetwork
            if numel(this.ActionInfo) > 1
                error(message('rl:agent:errSACMultiActionChannel'))
            end
            
            % inputs parsing
            switch nargin
                case 2
                    Option = rlSACAgentOptions;
                    Critic2 = [];
                case 3
                    validateattributes(varargin{1}, {'rl.representation.rlQValueRepresentation','rl.option.rlSACAgentOptions'}, {'scalar', 'nonempty'}, mfilename, '', 2);
                    if isa(varargin{1},'rl.option.rlSACAgentOptions')
                        Critic2 = [];
                        Option = varargin{1};
                    elseif isa(varargin{1},'rl.representation.rlQValueRepresentation')
                        Critic2 = varargin{1};
                        Option = rlSACAgentOptions;
                    end
                case 4
                    Critic2 = varargin{1};
                    Option = varargin{2};
                    validateattributes(Critic2, {'rl.representation.rlQValueRepresentation'}, {'scalar', 'nonempty'}, mfilename, 'Critic2', 3);
                    validateattributes(Option, {'rl.option.rlSACAgentOptions'}, {'scalar', 'nonempty'}, mfilename, 'Options', 4);
            end
            
            % set agent option
            this.AgentOptions = Option;
            
            % set representations
            setActor(this, Actor);
            setCritic(this, Critic1, Critic2);
            setEntropy(this);
            this.HasCritic = true;
            
        end
        function TargetEntropy = getTargetEntropy(this)
            % return the value of the TargetEntropy field
            TargetEntropy = exp(this.TargetLogEntropy);
        end
        function Actions = getActionWithExploration(this, Observation)
            Actions = getAction(this.Actor, Observation);
            Actions = this.ActionInfo.saturate(Actions);
            if numel(Actions) == 1
                Actions = Actions{1};
            end
        end
        function [critic1, critic2] = getCritic(this)
            % getCritic: Return the critic representation, CRITIC1 and
            % CRITIC2 for the SAC agent, SACAgent
            %
            % [CRITIC1, CRITIC2] = getCritic(SACAGent)
            critic1 = this.Critic1;
            if DoUseCritic2(this)
                critic2 = this.Critic2;
            else
                critic2 = [];
            end
        end
        function this = setCritic(this, Critic1, Critic2)
            % validate critic is a single output Q representation
            validateattributes(Critic1, {'rl.representation.rlQValueRepresentation'}, {'scalar', 'nonempty'}, '', 'Critic1');
            
            % validate critic is single output Q representations
            if strcmpi(getQType(Critic1),'multiOutput')
                error(message('rl:agent:errSACMultiQ'))
            end
            
            % validate action and observation infos are same
            checkAgentRepDataSpecCompatibility(this, Critic1)
            
            % set critic loss
            Critic1 = setLoss(Critic1,"rl.loss.rlmse");
            
            % Set critic network
            this.Critic1 = Critic1;
            
            % Construct target network
            this.TargetCritic1 = this.Critic1;
            
            if ~isempty(Critic2)
                validateattributes(Critic2, {'rl.representation.rlQValueRepresentation'}, {'scalar', 'nonempty'}, '', 'Critic2');
                if strcmpi(getQType(Critic2),'multiOutput')
                    error(message('rl:agent:errSACMultiQ'))
                end
                checkAgentRepDataSpecCompatibility(this, Critic2)
                
                if isequal(Critic1, Critic2)
                    % if 2 critics are actually 1 handle object, or critic 
                    % 2 is a deep copy of critic 1 with same params
                    error(message('rl:agent:errDoubleQSame'))
                end
                
                Critic2 = setLoss(Critic2,"rl.loss.rlmse");
                this.Critic2 = Critic2;
                this.TargetCritic2 = this.Critic2;
            end
            reset(this);
        end
        function Actor = getActor(this)
            % getActor: Return the actor representation, ACTOR, for the
            % specified reinforcement learning agent.
            %
            %   ACTOR = getActor(AGENT)
            
            Actor = this.Actor;
        end
        function this = setActor(this, Actor)
            % setActor: Set the actor of the reinforcement learning agent
            % using the specified representation, ACTOR, which must be
            % consistent with the observations and actions of the agent.
            %
            %   AGENT = setActor(AGENT,ACTOR)
            
            % validate actor is a deterministic actor representation
            validateattributes(Actor, {'rl.representation.rlStochasticActorRepresentation'}, {'scalar', 'nonempty'}, '', 'Actor');
            
            % validate action and observation infos are same
            checkAgentRepDataSpecCompatibility(this, Actor)
            
            % set SAC loss
            Actor = setLoss(Actor, "rl.loss.sacNoTanh");
            
            % Set actor network
            this.Actor = Actor;
            
            reset(this)
        end
        function set.AgentOptions(this,NewOptions)
            validateattributes(NewOptions,{'rl.option.rlSACAgentOptions'},{'scalar'},'','AgentOptions');
            
            % check if the experience buffer needs to be rebuild
            rebuildExperienceBuffer = isempty(this.ExperienceBuffer) || ...
                this.AgentOptions_.ExperienceBufferLength ~= NewOptions.ExperienceBufferLength;
            % rebuildEntropy = isempty(this.AgentOptions_) || ~isequal(this.AgentOptions_.EntropyWeightOptions, NewOptions.EntropyWeightOptions);
            this.AgentOptions_ = NewOptions;
            this.SampleTime = NewOptions.SampleTime;
            
            % build the experience buffer if necessary
            if rebuildExperienceBuffer
                if isempty(this.ExperienceBuffer)
                    buildBuffer(this);
                else
                    resize(this.ExperienceBuffer,this.AgentOptions_.ExperienceBufferLength);
                end
                this.ExperienceBuffer.DoValidate = false;
            end
            % save the experience buffer with the agent
            setSaveMemoryWithBuffer(this.ExperienceBuffer,this.AgentOptions_.SaveExperienceBufferWithAgent);
            
            % change entropy related properties
            %if rebuildEntropy
                setEntropy(this);
            %end
        end
        function Options = get.AgentOptions(this)
            Options = this.AgentOptions_;
        end
    end
    %======================================================================
    % Implementation of abstract methods
    %======================================================================
    methods (Access = protected)
        function [rep,argStruct] = generateProcessFunctions_(this,argStruct)
            % SAC can use the default code gen
            rep = this.Actor;
        end
        function Actions = getActionImpl(this, Observation)
            % Return mean values for the observations
            % REVIST: does not support multi-action channel
            ProbabilityParams = evaluate(this.Actor, Observation);
            ProbabilityParams = ProbabilityParams{1};
            NumActions = numel(this.ActionInfo);
            Actions = ProbabilityParams(1:NumActions, :);
            Actions = this.ActionInfo.saturate(Actions);
            if numel(Actions) == 1
                Actions = Actions{1};
            end
        end
        function resetImpl(this)
            % rebuild agent properties due to any potential changes in
            % options
            if this.AgentOptions.ResetExperienceBufferBeforeTraining        
                % Construct replay memory
                buildBuffer(this);
                this.ExperienceBuffer.DoValidate = false;
            end
            this.TargetUpdateCounter = 1;
        end
        function Q0 = evaluateQ0(this, exp)
            Observation = exp{1};
            Action = getAction(this.Actor, Observation);
            Action = this.ActionInfo.saturate(Action);
            Q0 = getValue(this.Critic1, Observation, Action);
            if DoUseCritic2(this)
                Q0 = min( Q0, getValue(this.Critic2, Observation, Action));
            end
            if isa(Q0,'dlarray')
                Q0 = extractdata(Q0);
            end
        end
        %set/get tunable parameters
        function setLearnableParametersImpl(this, p)
            this.Actor = setLearnableParameters(this.Actor, p.Actor);
            this.Critic1 = setLearnableParameters(this.Critic1,p.Critic1);
            if DoUseCritic2(this)
                this.Critic2 = setLearnableParameters(this.Critic2, p.Critic2);
            end
        end
        
        function p = getLearnableParametersImpl(this)
            p.Actor  = getLearnableParameters(this.Actor);
            p.Critic1 = getLearnableParameters(this.Criti1c);
            if DoUseCritic2(this)
                p.Critic2 = getLearnableParameters(this.Critic2);
            end
        end
        function Action = learn(this,exp)
            % learn from the current set of experiences where
            % exp = {state,action,reward,nextstate,isdone}
            % Return the noisy action.
            % NOTE: only learn, update params and exploration once exp
            % buffer has more experiences than MiniBatchSize
            
            % store experiences
            appendExperience(this,exp);
            
            % generate a minibatch: MiniBatchSize length of cell array with
            % {state,action,reward,nextstate,isdone} elements
            minibatch = createSampledExperienceMiniBatch(...
                this.ExperienceBuffer,...
                this.AgentOptions.MiniBatchSize,...
                this.AgentOptions.DiscountFactor,...
                this.AgentOptions.NumStepsToLookAhead);
                        
            if ~isempty(minibatch)
                % perform the learning on the representation
                stepRepresentation(this,minibatch)
            end
            
            % compute action from the current policy
            % {state,action,reward,nextstate,isdone}
            Action = getActionWithExploration(this,exp{4});
        end
        function trainingOptions = validateAgentTrainingCompatibilityImpl(this,trainingOptions)
            % Validate SAC agent training options compatibility
            if ~strcmpi(trainingOptions.Parallelization,'none')
                dataToSend = trainingOptions.ParallelizationOptions.DataToSendFromWorkers;
                % SAC agent only support send experiences for parallel
                if ~strcmpi(dataToSend,'Experiences')
                    error(message('rl:general:errParallelSendGradNotSupport'));
                end
                
                % currently only support NumStepsToLookAhead = 1
                if this.AgentOptions.NumStepsToLookAhead ~= 1
                    error(message('rl:general:errParallelSendExpNumStepLookAheadEQ1'));
                end
            end
        end
        function HasState = hasStateImpl(this)
            HasState = hasState(this.Actor);
        end
        function that = copyElement(this)
            that = copyElement@rl.agent.AbstractAgent(this);
            that.ExperienceBuffer = copy(this.ExperienceBuffer);
        end
    end
    
    methods (Hidden)
        function preSettings = preTrain(this)
            preSettings = preTrain@rl.agent.AbstractAgent(this);
            preSettings.ValidateInputArgumentsForExperienceBuffer = ...
                this.ExperienceBuffer.DoValidate;
            this.ExperienceBuffer.DoValidate = false;
        end
        
        function postTrain(this,preSettings)
            postTrain@rl.agent.AbstractAgent(this,preSettings);
            this.ExperienceBuffer.DoValidate = ...
                preSettings.ValidateInputArgumentsForExperienceBuffer;
        end
        
        function appendExperience(this,experiences)
            % append experiences to buffer
            append(this.ExperienceBuffer,{experiences});
        end
        
        function applyGradient(this,g)
            % update representation from gradients
            if ~isempty(g)
                this.Critic1 = optimize(this.Critic1,g.Critic1);
                if DoUseCritic2(this)
                    this.Critic2 = optimize(this.Critic2,g.Critic2);
                end
                if ~isempty(g.Actor)
                    this.Actor = optimize(this.Actor ,g.Actor);
                end
                updateTargetRepresentations(this);
            end
        end
        
        % Accumulate representation gradients from current experiences,
        % return gradient values and parameter values
        function g = accumulateGradient(this,~)
            % generate a minibatch: MiniBatchSize length of cell array with
            % {state,action,reward,nextstate,isdone} elements
            minibatch = createSampledExperienceMiniBatch(...
                this.ExperienceBuffer,...
                this.AgentOptions.MiniBatchSize,...
                this.AgentOptions.DiscountFactor,...
                this.AgentOptions.NumStepsToLookAhead);
                        
            if isempty(minibatch)
                g = [];
            else
                % perform the learning on the representation
                g = stepRepresentation(this,minibatch);
            end
        end
    end
    
    %======================================================================
    % Step representation methods
    %======================================================================
    methods (Access = private)
        function updateTargetRepresentations(this)
            % Update the target networks
            this.TargetCritic1 = updateTarget(this,this.Critic1,this.TargetCritic1,...
                this.AgentOptions.TargetSmoothFactor,...
                this.AgentOptions.TargetUpdateFrequency);
            if DoUseCritic2(this)
                this.TargetCritic2 = updateTarget(this,this.Critic2,this.TargetCritic2,...
                this.AgentOptions.TargetSmoothFactor,...
                this.AgentOptions.TargetUpdateFrequency);
            end
        end
        function this = setEntropy(this)
            if this.AgentOptions.EntropyWeightOptions.TargetEntropy == -1
                this.TargetLogEntropy = -numel(this.ActionInfo);
            else
                this.TargetLogEntropy = log(this.AgentOptions.EntropyWeightOptions.TargetEntropy);
            end
            this.EntropyWeight = dlarray(this.AgentOptions.EntropyWeightOptions.EntropyWeight);
            this.LogEntropyWeight = log(this.EntropyWeight);
            this.LearnableEntropy = true;
            if this.AgentOptions.EntropyWeightOptions.LearnRate == 0
                this.LearnableEntropy = false;
            end
            if this.LearnableEntropy
                EntropyOptimizerName = this.AgentOptions.EntropyWeightOptions.Optimizer;
                EntropyOptimizerOptions = this.AgentOptions.EntropyWeightOptions;
                
                switch EntropyOptimizerName
                    case "adam"
                        this.EntropyOptimizer = rl.internal.optimizer.rlADAMSolver(EntropyOptimizerOptions);
                    case "sgdm"
                        this.EntropyOptimizer = rl.internal.optimizer.rlSGDMSolver(EntropyOptimizerOptions);
                    case "rmsprop"
                        this.EntropyOptimizer = rl.internal.optimizer.rlRMSPropSolver(EntropyOptimizerOptions);
                    otherwise
                        error(message('rl:agent:errSACEntropyOptim'));
                end
            end
        end
        function this = setActionScale(this)
            LowerLimit = this.ActionInfo.LowerLimit;
            UpperLimit = this.ActionInfo.UpperLimit;
            Scale = (UpperLimit - LowerLimit) ./ 2;
            Shift = UpperLimit - Scale;
            this.ActionScale = Scale;
            this.ActionShift = Shift;
        end
        function varargout = stepRepresentation(this,minibatch)
            % Update the critic, actor, target critic and target actor
            % given a minibatch and a target discount factor
            
            if nargout
                % REVISIT Gorilla parallel training for SAC (send gradients for off-policy agent)
                
                % This branch supports gorilla parallel training where
                % workers send gradients to learners. 'train' methods
                % output a struct that contains critic and actor gradients.
                % The gradient will be supply to applyGradient(this,grad)
                [s.Critic, s.Critic2] = trainCriticWithBatch(this,minibatch);
                s.Actor  = trainActorWithBatch(this,minibatch);
                varargout{1} = s;
            else
                % update the critic
                trainCriticWithBatch(this,minibatch);
                % update the actor
                trainActorWithBatch(this,minibatch);
                % update entropy weight
                trainEntropyWeight(this,minibatch);
                % update the target representations
                updateTargetRepresentations(this);
                
            end
        end
        function varargout = trainCriticWithBatch(this, miniBatch)
            % update the critic against a minibatch set
            
            Observations     = miniBatch{1};
            Actions          = miniBatch{2};
            Rewards          = miniBatch{3};
            NextObservations = miniBatch{4};
            IsDones          = miniBatch{5};
           
            ProbabilityParams = evaluate(this.Actor, NextObservations);
            NextActions = sample(this.Actor.SamplingStrategy, ProbabilityParams);
            AdjustedNextActions = this.ActionInfo.saturate(NextActions);
            TargetQ = getValue(this.TargetCritic1, NextObservations, AdjustedNextActions);
            if DoUseCritic2(this)
                TargetQ2 = getValue(this.TargetCritic2, NextObservations, AdjustedNextActions);
                TargetQ = min(TargetQ, TargetQ2);
            end
            Densities = evaluate(this.Actor.SamplingStrategy, ProbabilityParams{1}, AdjustedNextActions{1});
            Densities = prod(Densities, 1);
            LogDensities = log(rl.internal.dataTransformation.boundAwayFromZero(Densities));
            TargetQ = TargetQ - this.EntropyWeight * LogDensities;
            
            DoneIdx = IsDones == 1;
            Gamma = this.AgentOptions.DiscountFactor;
            n = this.AgentOptions.NumStepsToLookAhead;
            
            % get target Q values we should expect the network to work
            % towards
            TargetQ(~DoneIdx) = Rewards(~DoneIdx) + (Gamma^n).*TargetQ(~DoneIdx);
            
            % for final step, just use the immediate reward, since there is
            % no more a next state
            TargetQ(DoneIdx) = Rewards(DoneIdx);
            % train the critic or get the gradients
            CriticGradient1 = gradient(this.Critic1,'loss-parameters',...
                    [Observations, Actions], TargetQ);
            if DoUseCritic2(this)
                CriticGradient2 = gradient(this.Critic2,'loss-parameters',...
                    [Observations, Actions], TargetQ);
            end
            if nargout
                s.Critic1 = CriticGradient1;
                if DoUseCritic2(this)
                    s.Critic2 = CriticGradient2;
                else
                    s.Critic2 = [];
                end
            else
                this.Critic1 = optimize(this.Critic1, CriticGradient1);
                if DoUseCritic2(this)
                    this.Critic2 = optimize(this.Critic2, CriticGradient2);
                end
            end
            if nargout
                varargout{1} = s;
            end
        end
        function varargout = trainActorWithBatch(this, miniBatch)
            % update the actor against a minibatch set
            LossVariable.Observations = miniBatch{1};
            LossVariable.Actor = this.Actor;
            LossVariable.Critic1 = this.Critic1;
            LossVariable.Critic2 = this.Critic2;
            LossVariable.EntropyWeight = this.EntropyWeight;
            LossVariable.ActionInfo = this.ActionInfo;
            % pass observation and loss varaible to loss
            ActorGradient = gradient(this.Actor, 'loss-parameters', LossVariable.Observations, LossVariable);
            if nargout
                varargout{1} = ActorGradient;
            else
                this.Actor = optimize(this.Actor, ActorGradient);
            end
        end
        function trainEntropyWeight(this, miniBatch)
            if this.LearnableEntropy
                Observations = miniBatch{1};
                ProbabilityParams = evaluate(this.Actor, Observations);
                Actions = sample(this.Actor.SamplingStrategy, ProbabilityParams);
                AdjustedActions = this.ActionInfo.saturate(Actions);
                Densities = evaluate(this.Actor.SamplingStrategy, ProbabilityParams{1}, AdjustedActions{1});
                Densities = prod(Densities, 1);
                LogDensities = log(rl.internal.dataTransformation.boundAwayFromZero(Densities));
                EntropyGradient = - mean( LogDensities + this.TargetLogEntropy, 'all');
                EntropyGradient = max(EntropyGradient, -this.AgentOptions.EntropyWeightOptions.GradientThreshold);
                EntropyGradient = min(EntropyGradient, this.AgentOptions.EntropyWeightOptions.GradientThreshold);
                [this.EntropyOptimizer, this.LogEntropyWeight] = calculateUpdate(this.EntropyOptimizer, this.LogEntropyWeight, ...
                    EntropyGradient, this.AgentOptions.EntropyWeightOptions.LearnRate);
                this.EntropyWeight = exp(this.LogEntropyWeight);
            end
        end
    end
    methods (Access = private)
        % TODO: should these be an util function or protected for DDPG, TD3,
        % SAC subclass
        function buildBuffer(this)
            this.ExperienceBuffer = rl.util.ExperienceBuffer(...
                this.AgentOptions_.ExperienceBufferLength, ...
                this.ObservationInfo, ...
                this.ActionInfo);
        end
        function checkAgentRepDataSpecCompatibility(this, Rep)
            if ~isCompatible(this.ActionInfo, Rep.ActionInfo)
                error(message('rl:agent:errActionInfoAC'))
            end
            if ~isCompatible(this.ObservationInfo, Rep.ObservationInfo)
                error(message('rl:agent:errObservationInfoAC'))
            end
        end
    end
    
    methods (Access = private)
        % TODO: should these be an util function or protected for any agent
        % that uses 2 critics?
        function tf = DoUseCritic2(this)
            tf = ~isempty(this.Critic2);
        end
    end
end