classdef rlTD3Agent < rl.agent.AbstractAgent
    % rlTD3Agent: Implements twin delayed deep deterministic policy gradient agent
    
    % Copyright 2019 The MathWorks Inc.
    
    properties (Dependent)
        % Options to configure RL agent
        AgentOptions
    end
    
    properties (Access = private)
        % Private options to configure RL agent
        AgentOptions_ = [];
        % Counter determines when to update the policy
        PolicyUpdateCounter
    end
    
    properties (SetAccess = private)
        % Circular buffer for experience replay
        ExperienceBuffer
    end
    
    properties (Access = private)
        % Critic function approximator
        Critic
        
        % Target critic function approximator
        TargetCritic
        
        % Actor function approximator
        Actor
        
        % Target Actor function approximator
        TargetActor
        
        % Noise model
        ExplorationModel
        TargetPolicySmoothModel
    end
    
    methods
        function this = rlTD3Agent(Actor, Critic, Option)
            % Constructor
            
            this = this@rl.agent.AbstractAgent();
            
            % extract observation and action info
            this.ActionInfo = Actor.ActionInfo;
            this.ObservationInfo = Actor.ObservationInfo;
            
            % REVISIT: Not support multiple action channels due to DLT
            % limitation (a loss layer cannot take multiple inputs)
            % multi actions might work with dlnetwork
            if numel(this.ActionInfo) > 1
                error(message('rl:agent:errDPGMultiActionChannel'))
            end
            
            % set agent option
            this.AgentOptions = Option;
            
            % set representations
            setActor(this,Actor);
            setCritic(this,Critic);
            this.HasCritic = true;
        end
        
        % Action methods
        function action = getActionWithExploration(this,observation)
            % Given the current state of the system, return an action.
            action = getAction(this.Actor,observation);
            % perturb the action
            action = applyNoise(this.ExplorationModel,action);
            % saturate the actions
            action = saturate(this.ActionInfo, action);
            
            % return numeric if action is single channel for backward
            % compatibility
            if numel(action) == 1
                action = action{1};
            end
        end
        
        %==================================================================
        % Get/set
        %==================================================================
        function Critic = getCritic(this)
            % getCritic: Return the critic representations vector, for the 
            % TD3 agent, TD3AGENT.
            %
            %   CRITIC = getCritic(TD3AGENT)
            
            Critic = this.Critic;
        end
        function this = setCritic(this,Critic)
            % setCritic: Set the critic of the reinforcement learning agent
            % using the specified representation, CRITIC, which must be
            % consistent with the observations and actions of the agent.
            %
            %   AGENT = setCritic(AGENT,CRITIC)
            
            % validate critic is a single output Q representation
            validateattributes(Critic, {'rl.representation.rlQValueRepresentation'}, {'vector', 'nonempty'}, '', 'Critic');
            
            if numel(Critic) > 2
                error(message('rl:agent:errTD3NumCriticGt2'))
            end
            
            for ct = 1:numel(Critic)
                % REVISIT: not support agent has state
                if hasState(Critic(ct))
                    error(message('rl:agent:errAgentHasStateNotSupport'))
                end
                
                % validate critics are single output Q representations
                if strcmpi(getQType(Critic(ct)),'multiOutput')
                    error(message('rl:agent:errDPGMultiQ'))
                end
                
                % validate action and observation infos are same
                checkAgentRepDataSpecCompatibility(this, Critic(ct));
                
                % set critic loss
                Critic(ct) = setLoss(Critic(ct),"rl.loss.rlmse");
                
                if ct > 1 && isequal(getLearnableParameters(Critic(ct)), getLearnableParameters(Critic(1))) ...
                        && isequal(Critic(ct), Critic(1))
                    % REVISIT: if numel(Critic)>2, need to check mutual
                    error(message('rl:agent:errDoubleQSame'))
                end
            end
            
            % set critic
            this.Critic = Critic;
                
            % construct target network
            this.TargetCritic = Critic;
            
            reset(this);
        end
        function actor = getActor(this)
            % getActor: Return the actor representation, ACTOR, for the
            % specified reinforcement learning agent.
            %
            %   ACTOR = getActor(AGENT)
            %
            
            actor = this.Actor;
        end
        function this = setActor(this, Actor)
            % setActor: Set the actor of the reinforcement learning agent
            % using the specified representation, ACTOR, which must be
            % consistent with the observations and actions of the agent.
            %
            %   AGENT = setActor(AGENT,ACTOR)
            
            % validate actor is a deterministic actor representation
            validateattributes(Actor, {'rl.representation.rlDeterministicActorRepresentation'}, {'scalar', 'nonempty'}, '', 'Actor');
            
            % REVISIT: not support agent has state
            if hasState(Actor)
                error(message('rl:agent:errAgentHasStateNotSupport'))
            end
            
            % validate action and observation infos are same
            checkAgentRepDataSpecCompatibility(this, Actor)
            
            % set dpg loss
            Actor = setLoss(Actor,"rl.loss.dpg");
            
            % set actor network
            this.Actor = Actor;
            
            % construct target network
            this.TargetActor  = this.Actor;
            
            reset(this)
        end
        function set.AgentOptions(this,NewOptions)
            validateattributes(NewOptions,{'rl.option.rlTD3AgentOptions'},{'scalar'},'','AgentOptions');
            % TODO: switch between OU and Gaussian noise
            
            % check if the experience buffer needs to be rebuild
            rebuildExperienceBuffer = isempty(this.ExperienceBuffer) || ...
                this.AgentOptions_.ExperienceBufferLength ~= NewOptions.ExperienceBufferLength;
            % check to see if we need to rebuild the noise model
            rebuildNoise = isempty(this.ExplorationModel) || ...
                ~isequal(this.AgentOptions_.ExplorationModel,NewOptions.ExplorationModel);
            
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
            
            % build the noise model if necessary
            if rebuildNoise
                % extract the noise options
                explorationNoiseOpts = this.AgentOptions_.ExplorationModel;
                targetPolicyNoiseOpts = this.AgentOptions_.TargetPolicySmoothModel;
                % create the noise model
                actionDims = {this.ActionInfo.Dimension}';
                this.ExplorationModel = rl.util.createNoiseModelFactory(...
                    actionDims,explorationNoiseOpts,getSampleTime(this));
                this.TargetPolicySmoothModel = rl.util.createNoiseModelFactory(...
                    actionDims,targetPolicyNoiseOpts,getSampleTime(this));
            end
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
            % TD3 can use the default code gen
            
            rep = this.Actor;
        end
        
        function action = getActionImpl(this, observation)
            % Given the current state of the system, return a saturated
            % action.
            
            action = getAction(this.Actor, observation);
            action = saturate(this.ActionInfo, action);
            
            % return numeric if action is single channel for backward
            % compatibility
            if numel(action) == 1
                action = action{1};
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
            
            % reset the noise models
            reset(this.ExplorationModel);
            reset(this.TargetPolicySmoothModel);
            
            % utils to keep track of update
            this.PolicyUpdateCounter = 1;
        end
        
        function q0 = evaluateQ0(this,exp)
            % overload for agents that implement critics
            observation = exp{1};
            action = getAction(this, observation);
            for ct = 1:numel(this.Critic)
                q0(ct) = getValue(this.Critic(ct), observation ,action); %#ok<AGROW>
            end
            q0 = min(q0);
            if isa(q0,'dlarray')
                q0 = extractdata(q0);
            end
        end
        
        % set/get tunable parameters
        function setLearnableParametersImpl(this,p)
            this.Actor = setLearnableParameters(this.Actor, p.Actor);
            for ct = 1:numel(this.Critic)
                this.Critic(ct) = setLearnableParameters(this.Critic(ct),p.Critic{ct});
            end
        end
        function p = getLearnableParametersImpl(this)
            p.Actor  = getLearnableParameters(this.Actor );
            for ct = 1:numel(this.Critic)
                p.Critic{ct} = getLearnableParameters(this.Critic(ct));
            end
        end
        
        function action = learn(this,exp)
            % learn from the current set of experiences where
            % exp = {state,action,reward,nextstate,isdone}
            % Return the noisy action.
            
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
            action = getActionWithExploration(this,exp{4});
        end
        
        function trainingOptions = validateAgentTrainingCompatibilityImpl(this,trainingOptions)
            % Validate TD3 agent training options compatibility
            
            if ~strcmpi(trainingOptions.Parallelization,'none')
                dataToSend = trainingOptions.ParallelizationOptions.DataToSendFromWorkers;
                % TD3 agent only support send experiences for parallel
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
            that.ExplorationModel = copy(this.ExplorationModel);
            that.TargetPolicySmoothModel = copy(this.TargetPolicySmoothModel);
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
                for ct = 1:numel(this.Critic)
                    this.Critic(ct) = optimize(this.Critic(ct),g.Critic{ct});
                end
                if ~isempty(g.Actor)
                    this.Actor = optimize(this.Actor ,g.Actor);
                end
                % REVISIT for parallel
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
            % Update the target representations
            
            for ct = 1:numel(this.Critic)
                this.TargetCritic(ct) = updateTarget(this,this.Critic(ct),this.TargetCritic(ct),...
                    this.AgentOptions.TargetSmoothFactor,...
                    this.AgentOptions.TargetUpdateFrequency);
            end
            this.TargetActor = updateTarget(this,this.Actor,this.TargetActor,...
                this.AgentOptions.TargetSmoothFactor,...
                this.AgentOptions.TargetUpdateFrequency);
        end
        
        % Step representation
        function varargout = stepRepresentation(this,minibatch)
            % update the critic, actor, target critic and target actor
            % given a minibatch and a target discount factor
            
            if nargout
                % get the gradients
                s.Critic = trainCriticWithBatch(this,minibatch);
                if mod(this.PolicyUpdateCounter, this.AgentOptions.PolicyUpdateFrequency) == 0
                    s.Actor = trainActorWithBatch(this,minibatch);
                else
                    this.PolicyUpdateCounter = this.PolicyUpdateCounter + 1;
                    s.Actor = [];
                end
                varargout{1} = s;
            else
                % update the actor and critic
                trainCriticWithBatch(this,minibatch);
                if mod(this.PolicyUpdateCounter, this.AgentOptions.PolicyUpdateFrequency) == 0
                    trainActorWithBatch(this,minibatch);
                else
                    this.PolicyUpdateCounter = this.PolicyUpdateCounter + 1;
                end
                % update the target networks
                updateTargetRepresentations(this);
            end
        end
        
        function varargout = trainCriticWithBatch(this,miniBatch)
            % update the critics and actor against a minibatch set
            
            observations     = miniBatch{1};
            actions          = miniBatch{2};
            rewards          = miniBatch{3};
            nextobservations = miniBatch{4};
            isdones          = miniBatch{5};
            
            % compute the next actions from the target actor
            nextactions = getAction(this.TargetActor,nextobservations);
            % perturb the action with target noise
            nextactions = applyNoise(this.TargetPolicySmoothModel,nextactions);
            % saturate the actions
            nextactions = saturate(this.ActionInfo, nextactions);
            
            % compute the next step expected Q value (bootstrapping)
            for ct = 1:numel(this.Critic)
                if ct < 2
                    targetq = getValue(this.TargetCritic(ct), nextobservations, nextactions);
                else
                    targetq = min(targetq, getValue(this.TargetCritic(ct), nextobservations, nextactions));
                end
            end
            
            doneidx = isdones == 1;
            gamma = this.AgentOptions.DiscountFactor;
            n = this.AgentOptions.NumStepsToLookAhead;
            
            % get target Q values we should expect the network to work
            % towards
            targetq(~doneidx) = rewards(~doneidx) + (gamma^n).*targetq(~doneidx);
            
            % for final step, just use the immediate reward, since there is
            % no more a next state
            targetq(doneidx) = rewards(doneidx);
            
            % train the critic or get the gradients
            for ct = 1:numel(this.Critic)
                criticGradient(ct) = {gradient(this.Critic(ct),'loss-parameters',...
                    [observations, actions], targetq)}; %#ok<AGROW>
            end
            if nargout
                s.Critic = criticGradient;
            else
                for ct = 1:numel(this.Critic)
                    this.Critic(ct) = optimize(this.Critic(ct), criticGradient{ct});
                end
            end
            if nargout
                varargout{1} = s;
            end
        end
        
        function varargout = trainActorWithBatch(this,miniBatch)
            % update the actor against a minibatch set
            
            
            lossVariable.Observation = miniBatch{1};
            lossVariable.Critic = this.Critic;
            
            % pass observation and critic to 'rl.loss.dpg' loss
            actorGradient = gradient(this.Actor,'loss-parameters',lossVariable.Observation,lossVariable);
            if nargout
                % REVISIT for parallel training
                varargout{1} = actorGradient;
            else
                this.Actor = optimize(this.Actor, actorGradient);
            end
            
            % reset policy update counter
            this.PolicyUpdateCounter = 1;
        end
    end
    
    methods (Access = private)
        % REVISIT: should these be an util function or protected for DDPG, TD3,
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
end