classdef rlDDPGAgent < rl.agent.AbstractAgent
    % rlDDPGAgent: Implements deep deterministic policy gradient agent
    % TODO: error checking support only single action channel
    
    % Copyright 2017-2018 The MathWorks Inc.

    properties (Dependent)
        % Options to configure RL agent
        AgentOptions
    end
    
    properties (Access = private)
        % Private options to configure RL agent
        AgentOptions_ = [];
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
        NoiseModel
    end
    
    properties (Constant, Access = private)
        % version indicator for backward compatibility
        Version = 2
    end

    methods
        function this = rlDDPGAgent(Actor, Critic, Option)
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
            setActor(this, Actor);
            setCritic(this, Critic);
            this.HasCritic = true;
        end
        
        % Action methods
        function action = getActionWithExploration(this,observation)
            % Given the current state of the system, return an action.
            
            action = getAction(this.Actor,observation);
            % perturb the action 
            action = applyNoise(this.NoiseModel,action);
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
        function critic = getCritic(this)
            % getCritic: Return the critic representation, CRITIC, for the
            % specified reinforcement learning agent, AGENT.
            %
            %   CRITIC = getCritic(AGENT)
            
            critic = this.Critic;
        end
        function this = setCritic(this, Critic)
            % setCritic: Set the critic of the reinforcement learning agent
            % using the specified representation, CRITIC, which must be
            % consistent with the observations and actions of the agent.
            %
            %   AGENT = setCritic(AGENT,CRITIC)
            
            % validate critic is a single output Q representation
            validateattributes(Critic, {'rl.representation.rlQValueRepresentation'}, {'scalar', 'nonempty'}, '', 'Critic');
            
            % REVISIT: not support agent has state
            if hasState(Critic)
                error(message('rl:agent:errAgentHasStateNotSupport'))
            end
            
            % validate critic is single output Q representations
            if strcmpi(getQType(Critic),'multiOutput')
                error(message('rl:agent:errDPGMultiQ'))
            end
            
            % validate action and observation infos are same
            checkAgentRepDataSpecCompatibility(this, Critic)
            
            % set critic loss
            Critic = setLoss(Critic,"rl.loss.rlmse");
            
            % Set critic network
            this.Critic = Critic;
            
            % Construct target network
            this.TargetCritic = this.Critic;
            
            reset(this);
        end
        function actor = getActor(this)
            % getActor: Return the actor representation, ACTOR, for the
            % specified reinforcement learning agent.
            %
            %   ACTOR = getActor(AGENT)
            
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
            
            % Set actor network
            this.Actor = Actor;
            
            % Construct target network
            this.TargetActor  = this.Actor;
            
            reset(this)
        end        
        function set.AgentOptions(this,NewOptions)
            validateattributes(NewOptions,{'rl.option.rlDDPGAgentOptions'},{'scalar'},'','AgentOptions');
            
            % check if the experience buffer needs to be rebuild
            rebuildExperienceBuffer = isempty(this.ExperienceBuffer) || ...
                this.AgentOptions_.ExperienceBufferLength ~= NewOptions.ExperienceBufferLength;
            % check to see if we need to rebuild the noise model
            rebuildNoise = isempty(this.NoiseModel) || ...
                ~isequal(this.AgentOptions_.NoiseOptions,NewOptions.NoiseOptions);
            
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
                noiseOpts = this.AgentOptions_.NoiseOptions;

                % create the noise model
                actionDims = {this.ActionInfo.Dimension}';
                this.NoiseModel = rl.util.createNoiseModelFactory(...
                    actionDims,noiseOpts,getSampleTime(this));
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
            % DDPG can use the default code gen
            
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
            
            % reset the noise model
            reset(this.NoiseModel);
        end
        
        function q0 = evaluateQ0(this,exp)
            % overload for agents that implement critics
            
            observation = exp{1};
            action = getAction(this, observation);
            q0 = getValue(this.Critic, observation, action);
            if isa(q0,'dlarray')
                q0 = extractdata(q0);
            end
        end
        
        % set/get tunable parameters
        function setLearnableParametersImpl(this,p)
            this.Actor = setLearnableParameters(this.Actor ,p.Actor );
            this.Critic = setLearnableParameters(this.Critic,p.Critic);
        end
        
        function p = getLearnableParametersImpl(this)
            p.Actor  = getLearnableParameters(this.Actor );
            p.Critic = getLearnableParameters(this.Critic);
        end
        
        function action = learn(this,exp)
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
            action = getActionWithExploration(this,exp{4});
        end
        
        function trainingOptions = validateAgentTrainingCompatibilityImpl(this,trainingOptions)
            % Validate DDPG agent training options compatibility
            
            if ~strcmpi(trainingOptions.Parallelization,'none')
                dataToSend = trainingOptions.ParallelizationOptions.DataToSendFromWorkers;
                % DDPG agent only support send experiences for parallel
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
            that.NoiseModel = copy(this.NoiseModel);
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
                this.Critic = optimize(this.Critic,g.Critic);
                this.Actor = optimize(this.Actor ,g.Actor);
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
            this.TargetCritic = updateTarget(this,this.Critic,this.TargetCritic,...
                this.AgentOptions.TargetSmoothFactor,...
                this.AgentOptions.TargetUpdateFrequency);
            this.TargetActor = updateTarget(this,this.Actor,this.TargetActor,...
                this.AgentOptions.TargetSmoothFactor,...
                this.AgentOptions.TargetUpdateFrequency);
        end
        
        function varargout = stepRepresentation(this,minibatch)
            % Update the critic, actor, target critic and target actor
            % given a minibatch and a target discount factor
            
            if nargout
                % This branch supports gorilla parallel training where
                % workers send gradients to learners. 'train' methods
                % output a struct that contains critic and actor gradients.
                % The gradient will be supply to applyGradient(this,grad)
                s.Critic = trainCriticWithBatch(this,minibatch);
                s.Actor  = trainActorWithBatch (this,minibatch);
                varargout{1} = s;
            else
                % update the critic
                trainCriticWithBatch(this,minibatch);
                % update the actor
                trainActorWithBatch(this,minibatch);
                % update the target representations
                updateTargetRepresentations(this);
            end
        end
        
        function varargout = trainCriticWithBatch(this,miniBatch)
            % update the critic against a minibatch set
            
            observations     = miniBatch{1};
            actions          = miniBatch{2};
            rewards          = miniBatch{3};
            nextobservations = miniBatch{4};
            isdones          = miniBatch{5};
            
            % compute the next actions from the target actor
            nextactions = getAction(this.TargetActor,nextobservations);
            
            % compute the next step expected Q value (bootstrapping)
            targetq = getValue(this.TargetCritic,nextobservations,nextactions);
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
            criticGradient = gradient(this.Critic,'loss-parameters',...
                    [observations, actions], targetq);
            if nargout
                varargout{1} = criticGradient;
            else
                this.Critic = optimize(this.Critic, criticGradient);
            end
        end
        % Actor methods
        function varargout = trainActorWithBatch(this,miniBatch)
            % update the actor against a minibatch set
            
            lossVariable.Observation = miniBatch{1};
            lossVariable.Critic = this.Critic;
            % pass observation and critic to 'rl.loss.dpg' loss
            actorGradient = gradient(this.Actor,'loss-parameters',lossVariable.Observation,lossVariable);
            if nargout
                varargout{1} = actorGradient;
            else
                this.Actor = optimize(this.Actor, actorGradient);
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
    
    methods (Static)
        function obj = loadobj(s)
            if isstruct(s)
                if ~isfield(s,'Version')
                    % version 1 does not have Version field
                    % In version 2,
                    %   - Critic changes from rlRepresentation to rlQValueRepresentation
                    %   - Actor changes from rlRepresentation to rlDeterministicActorRepresentation
                    %   - ExperienceBuffer requires obs and act dims inputs
                    criticModel = s.Critic.getModel;
                    targetCriticModel = s.TargetCritic.getModel;
                    if isa(criticModel,'DAGNetwork') || isa(criticModel,'nnet.cnn.LayerGraph')
                        s.Critic = rlQValueRepresentation(criticModel,...
                            s.ObservationInfo, s.ActionInfo, ...
                            'Observation', s.Critic.ObservationNames, ...
                            'Action', s.Critic.ActionNames, ...
                            s.Critic.Options);
                        s.TargetCritic = rlQValueRepresentation(targetCriticModel,...
                            s.ObservationInfo, s.ActionInfo, ...
                            'Observation', s.TargetCritic.ObservationNames, ...
                            'Action', s.TargetCritic.ActionNames, ...
                            s.TargetCritic.Options);
                    else
                        s.Critic = rlQValueRepresentation(criticModel,...
                            s.ObservationInfo, s.ActionInfo, ...
                            s.Critic.Options);
                        s.TargetCritic = rlQValueRepresentation(targetCriticModel,...
                            s.ObservationInfo, s.ActionInfo, ...
                            s.TargetCritic.Options);
                    end
                    
                    actorModel = s.Actor.getModel;
                    targetActorModel = s.TargetActor.getModel;
                    if isa(actorModel,'DAGNetwork') || isa(actorModel,'nnet.cnn.LayerGraph')
                        s.Actor = rlDeterministicActorRepresentation(actorModel,...
                            s.ObservationInfo, s.ActionInfo, ...
                            'Observation', s.Actor.ObservationNames, ...
                            'Action', s.Actor.ActionNames, ...
                            s.Actor.Options);
                        s.TargetActor = rlDeterministicActorRepresentation(targetActorModel,...
                            s.ObservationInfo, s.ActionInfo, ...
                            'Observation', s.TargetActor.ObservationNames, ...
                            'Action', s.TargetActor.ActionNames, ...
                            s.TargetActor.Options);
                    else
                        s.Actor = rlDeterministicActorRepresentation(actorModel,...
                            s.ObservationInfo, s.ActionInfo, ...
                            s.Actor.Options);
                        s.TargetActor = rlDeterministicActorRepresentation(targetActorModel,...
                            s.ObservationInfo, s.ActionInfo, ...
                            s.TargetActor.Options);
                    end
                end
                obj = rl.agent.rlDDPGAgent(s.Actor,s.Critic,s.AgentOptions_);
                % REVISIT: how to retain the state of optimizer
                obj.NoiseModel = s.NoiseModel;
                obj.TargetActor = s.TargetActor;
                obj.TargetCritic = s.TargetCritic;
                
                if obj.AgentOptions_.SaveExperienceBufferWithAgent
                    % only load the experience buffer if
                    % SaveExperienceBufferWithAgent is true
                    obj.ExperienceBuffer = s.ExperienceBuffer;
                end
            else
                obj = s;
            end
        end
    end
end