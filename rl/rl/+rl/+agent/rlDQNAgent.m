classdef rlDQNAgent < rl.agent.AbstractAgent
    % rlDQNAgent: Implements deep Q Network agent
    
    % Copyright 2017-2019 The MathWorks Inc.        
    
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

        % Epsilon-Greedy Exploration parameters
        ExplorationModel
    end
    
    properties (Constant, Access = private)
        % version indicator for backward compatibility
        Version = 2
    end

    methods
        function this = rlDQNAgent(Critic, Option)
            % Constructor
            
            this = this@rl.agent.AbstractAgent();
            
            % extract observation and action info
            this.ActionInfo = Critic.ActionInfo;
            this.ObservationInfo = Critic.ObservationInfo;
            
            % set agent option
            this.AgentOptions = Option;

            % set representation
            setCritic(this,Critic);
            this.HasCritic = true;
            
            % Construct exploration model
            this.ExplorationModel = this.AgentOptions.EpsilonGreedyExploration;
        end
        
        function Action = getActionWithExploration(this,Observation) 
            % Return an action with exploration
            
            if rand < this.ExplorationModel.Epsilon
                Action = usample(this.ActionInfo);  
                if hasState(this.Critic)
                    % To update a hidden state of LSTM, we first get the
                    % current hidden state (output of getValue).
                    [~, State] = getValue(this.Critic, Observation);

                    % We update the hidden state using the current
                    % hidden state.
                    this.Critic = setState(this.Critic,State);
                end
            else
                % Hidden state is update in getActionImp
                Action = getAction(this,Observation); 
            end
            
            if iscell(Action) && numel(Action) == 1
                % unwrap cell to produce same behavior with previous
                % release
                Action = Action{1};
            end
           
            if strcmp(getStepMode(this),"sim-with-exploration")
                % (for parallel training) update the noise model on workers
                this.ExplorationModel = update(this.ExplorationModel);
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
                        
            % validate critic is a Q representation
            validateattributes(Critic, {'rl.representation.rlQValueRepresentation'}, {'scalar', 'nonempty'}, '', 'Critic');
            
            % validate agent options
            validateOptionImpl(this, this.AgentOptions, Critic)
            
            % validate action and observation infos are same
            if ~isCompatible(this.ActionInfo,Critic.ActionInfo)
                error(message('rl:agent:errIncompatibleActionInfo'))
            end
            if ~isCompatible(this.ObservationInfo,Critic.ObservationInfo)
                error(message('rl:agent:errIncompatibleObservationInfo'))
            end
            
            % set loss      
            if hasState(Critic)
                % This loss function is for DRQN (DQN that uses LSTM).
                Critic = setLoss(Critic,"rl.loss.drq");
            else
                % This loss function is for DQN that does not use a recurrent
                % neural network layer (LSTM).
                Critic = setLoss(Critic,"rl.loss.dq");
            end
            
            % DQN does not support critic with continuous action data spec
            if isa(Critic.ActionInfo, 'rl.util.rlNumericSpec')
                error(message('rl:agent:errDQNContinuousCritic'))
            end
            
            % Set critic network
            this.Critic = Critic;
            
            % Construct target network
            this.TargetCritic = this.Critic;
            
            reset(this);
        end
        function set.AgentOptions(this,NewOptions)
            validateattributes(NewOptions,{'rl.option.rlDQNAgentOptions'},{'scalar'},'','AgentOptions');
            validateOptionImpl(this, NewOptions, this.Critic)
            % check if the experience buffer needs to be rebuild
            rebuildExperienceBuffer = isempty(this.ExperienceBuffer) || ...
                this.AgentOptions_.ExperienceBufferLength ~= NewOptions.ExperienceBufferLength;
            
            this.AgentOptions_ = NewOptions;
            this.SampleTime = NewOptions.SampleTime;
            this.ExplorationModel = this.AgentOptions.EpsilonGreedyExploration;
            
            % build the experience buffer if necessary
            if rebuildExperienceBuffer
                if isempty(this.ExperienceBuffer)
                    this.ExperienceBuffer = rl.util.ExperienceBuffer(this.AgentOptions_.ExperienceBufferLength, this.ObservationInfo, this.ActionInfo);
                else
                    resize(this.ExperienceBuffer,this.AgentOptions_.ExperienceBufferLength);
                end
                this.ExperienceBuffer.DoValidate = false;
            end
            % save the experience buffer with the agent
            setSaveMemoryWithBuffer(this.ExperienceBuffer,this.AgentOptions_.SaveExperienceBufferWithAgent);
        end
        function Options = get.AgentOptions(this)
            Options = this.AgentOptions_;
        end
    end
    
    %======================================================================
    % Implementation of abstract methods
    %======================================================================
    methods(Access = protected)
        function [rep,argStruct] = generateProcessFunctions_(this,argStruct)
            % DQN code gen
            rep = this.Critic;
            
            if strcmpi(getQType(this.Critic), 'multiOutput')
                argStruct = rl.codegen.generateDiscreteCriticMultiOutQFcn(argStruct,this.ActionInfo);
            else
                argStruct = rl.codegen.generateDiscreteCriticSingleOutQFcn(argStruct,this.ActionInfo);
            end
        end
        
        function q0 = evaluateQ0(this,exp)
            % overload for agents that implement critics
            observation = exp{1};
            if hasState(this.Critic)
                %LSTM resets its hidden state.
                this.Critic = resetState(this.Critic); 
            end
            q0 = getMaxQValue(this.Critic, observation);
            if isa(q0,'dlarray')
                q0 = extractdata(q0);
            end
        end
        
        % set/get tunable parameters
        function setLearnableParametersImpl(this,p)
            this.Critic = setLearnableParameters(this.Critic,p.Critic);
        end
        function p = getLearnableParametersImpl(this)
            p.Critic = getLearnableParameters(this.Critic);
        end
        
        function Action = getActionImpl(this,Observation) 
            % Given the current state of the system, return an action.

            % REVISIT: DRQN. Reset hidden state of LSTM here at the beginning
            % of each epock. When we use sim, it does not call learn().
            % Thus, we need to reset state here. This is a temporary solution.
            % We will need to have a reset funciton in agent so that each
            % agent can use its own reset function at the beginning of each
            % epoch. We may move the reset to stepImpl().
            if hasState(this.Critic)
                if this.EpisodeInfo.StepsTaken==1
                    this.Critic = resetState(this.Critic);
                end
                [~,ActionIdx, State] = getMaxQValue(this.Critic, Observation); 
                Action = getElementValue(this.ActionInfo,ActionIdx);
                this.Critic = setState(this.Critic, State);
            else
                [~,ActionIdx] = getMaxQValue(this.Critic, Observation);
                Action = getElementValue(this.ActionInfo,ActionIdx);
            end
            
        end
        
        function resetImpl(this)
            if this.AgentOptions.ResetExperienceBufferBeforeTraining
                % Construct replay memory
                this.ExperienceBuffer = rl.util.ExperienceBuffer(this.AgentOptions_.ExperienceBufferLength,this.ObservationInfo, this.ActionInfo);
                this.ExperienceBuffer.DoValidate = false;
            end
            
            % Revert exploration model to original parameters
            this.ExplorationModel = this.AgentOptions.EpsilonGreedyExploration;
            this.Critic = resetState(this.Critic);
        end
        
        function action = learn(this,exp)
            % learn from the current set of experiences where
            % exp = {state,action,reward,nextstate,isdone}
            % Return the action with exploration.
            % NOTE: only learn, update params and exploration once exp
            % buffer has more experiences than MiniBatchSize
                        
            
            % REVISIT: DRQN. Reset hidden state of lstm here at the beginning of each TRAINING epock.
            % This is a temporary solution. We will need to have a reset funciton in agent so that each
            % agent can use its own reset function at the beginning of each
            % epoch. We might move the resert to stepImpl().
            if this.EpisodeInfo.StepsTaken==1
                this.Critic = resetState(this.Critic);
            end
               
            % store experiences
            appendExperience(this,exp);
            
            % sample data from exp buffer and compute gradients
            Grad = accumulateGradient(this);
            
            % update critic params and update exploration
            applyGradient(this,Grad);
            
            % compute action from the current policy
            % {state,action,reward,nextstate,isdone}
            action = getActionWithExploration(this,exp{4});   
        end
        
        function trainingOptions = validateAgentTrainingCompatibilityImpl(this,trainingOptions)
            % Validate DQN agent training options compatibility
            
            if ~strcmpi(trainingOptions.Parallelization,'none')
                dataToSend = trainingOptions.ParallelizationOptions.DataToSendFromWorkers;
                % DQN agent only support send experiences for parallel
                if ~strcmpi(dataToSend,'Experiences')
                    error(message('rl:general:errParallelSendGradNotSupport'));
                end

                % currently only support NumStepsToLookAhead = 1
                if this.AgentOptions.NumStepsToLookAhead ~= 1
                    error(message('rl:general:errParallelSendExpNumStepLookAheadEQ1'));
                end
            end

            if hasState(this.Critic)
                if this.AgentOptions.SequenceLength > trainingOptions.MaxStepsPerEpisode
                    error(message('rl:agent:errSequenceLengthLongerThanMaxStepsPerEpisode'));
                end
                if this.AgentOptions.MiniBatchSize >= trainingOptions.MaxEpisodes
                    error(message('rl:agent:errMiniBatchSizeGtMaxEpisodes'));
                end
            end
        end
        
        function HasState = hasStateImpl(this)
            HasState = hasState(this.Critic);
        end
        
        function that = copyElement(this)
            that = copyElement@rl.agent.AbstractAgent(this);
            that.ExperienceBuffer = copy(this.ExperienceBuffer);
        end
    end
    
    %======================================================================
    % Step representation methods
    %======================================================================
    methods (Access = private)
        function CriticGradient = trainCriticWithBatch(this,MiniBatch, MaskIdx)
            % update the critic against a minibatch set
            % MaskIdx is for DRQN. MaskIdx indicates actual inputs (true)
            % and padded inputs (false). MaskIdx is empty for DQN.         
            % MaskIdx is a tensor of size [1 x MiniBatchSize x SequenceLength].
            
            if this.AgentOptions.UseDoubleDQN
                % DoubleDQN: r + DiscountFactor*Q[s',a' = argmax(qNetwork(s'))]
                LossVariable.TargetCritic = this.Critic;
            else
                % DQN:       r + DiscountFactor*Qtarget[s',a' = argmax(targetNetwork(s'))]
                LossVariable.TargetCritic = this.TargetCritic;
            end
            LossVariable.Action = MiniBatch{2};
            LossVariable.Reward = MiniBatch{3};
            LossVariable.NextObservation = MiniBatch{4};
            LossVariable.DoneIdx = MiniBatch{5} == 1;
            LossVariable.Discount = this.AgentOptions.DiscountFactor ^ ...
                this.AgentOptions.NumStepsToLookAhead;
            LossVariable.ActionInfo = this.ActionInfo;
                        
            if hasState(this.Critic)
                %LSTM uses a mask to ignore padded inputs.
                LossVariable.MaskIdx = MaskIdx;            
                CurrentState = getState(this.Critic); %This saves the current hidden state in the current environment.
                this.Critic = resetState(this.Critic); % We reset hidden state for the training.
            end
            
            % compute gradient
            if strcmpi(getQType(this.Critic),'singleOutput')
                CriticGradient = gradient(this.Critic,'loss-parameters',...
                    [MiniBatch{1},MiniBatch{2}], LossVariable);
            else
                CriticGradient = gradient(this.Critic,'loss-parameters',...
                    MiniBatch{1}, LossVariable);
            end
            if hasState(this.Critic)
                % After the trainig, we recover the state using the saved hidden state.
                % Otherwise, the cirtic network always uses the initialized
                % hidden state to choose the actions. This currentState is
                % not from the training.
                this.Critic = setState(this.Critic, CurrentState);
            end
        end
        
        function updateTargetRepresentations(this)
            % Update target critic parameters
            this.TargetCritic = updateTarget(this,this.Critic,this.TargetCritic,...
                this.AgentOptions.TargetSmoothFactor,...
                this.AgentOptions.TargetUpdateFrequency);
        end
        
        function validateOptionImpl(this,NewOptions,Critic)
            % Check compatibility of options and critic representation
            
            % DQN does not support actions that are matrix or tensor since
            % getElementIndicationMatrix in rlFiniteSetSpec does not
            % support them. It supports only scalar and vector.
            if iscell(this.ActionInfo.Elements)
                OneAction = this.ActionInfo.Elements{1};
            else
                OneAction = this.ActionInfo.Elements(1);
            end
            if ~isscalar(OneAction) && ~isvector(OneAction)
                error(message('rl:agent:errDQNNotSupportMatrixTensorAction'))
            end

            if ~isempty(Critic)
                if ~isempty(this.Critic)
                    if hasState(this.Critic) && ~hasState(Critic)
                        % If the original critic network architecture has state
                        % (e.g. LSTM, DRQN algorithm)), but new critic
                        % architecture does not have states (DQN algorithm), it
                        % is necessary to create a new agent with
                        % sequenceLength=1
                        error(message('rl:agent:errCriticChangedtoStateless'))

                    elseif hasState(this.Critic)==false && hasState(Critic)==true
                        % If the original critic network architecture does
                        % not have states(DQN algorithm), but new critic
                        % architecture has states (e.g. LSTM, DRQN
                        % algorithm), it is necessary to create a new agent
                        % with sequenceLength>1
                        error(message('rl:agent:errCriticChangedtoRNN'))
                    end
                end
                
                % n-Step Q learning is not supported for LSTM.
                if NewOptions.NumStepsToLookAhead > 1 && hasState(Critic)
                    error(message('rl:agent:errDQNwithLSTMNotSupportNStepQ'))
                end

                %  SequenceLength must be greater than one for LSTM
                if NewOptions.SequenceLength == 1 && hasState(Critic)
                    error(message('rl:agent:errDRQNSequenceLengthEq1'))
                end

                %  SequenceLength must be one for DQN
                if NewOptions.SequenceLength > 1 && ~hasState(Critic)
                    error(message('rl:agent:errStatelessDQNSequenceLengthGt1'))
                end
            end
        end
    end
    
    methods (Hidden)
        function Grad = accumulateGradient(this,~)
            % Accumulate representation gradients from current experiences,
            % return gradient values and parameter values
            
            % generate a minibatch: MiniBatchSize length of cell array with
            % {state,action,reward,nextstate,isdone} elements
                        
            MaskIdx = {}; %DQN does not use masking
            if hasState(this.Critic)
                [Minibatch, MaskIdx] = createSampledExperienceMiniBatchSequence(...
                    this.ExperienceBuffer,...
                    this.AgentOptions.MiniBatchSize,...
                    this.AgentOptions.SequenceLength);
            else
                Minibatch = createSampledExperienceMiniBatch(...
                    this.ExperienceBuffer,...
                    this.AgentOptions.MiniBatchSize,...
                    this.AgentOptions.DiscountFactor,...
                    this.AgentOptions.NumStepsToLookAhead);                
            end
                                   
            if isempty(Minibatch)
                Grad = [];
            else
                % perform the learning on the representation
                Grad.Critic = trainCriticWithBatch(this, Minibatch, MaskIdx);
                
                % update exploration
                this.ExplorationModel = update(this.ExplorationModel);
            end
        end
        
        function applyGradient(this,Grad)
            % Update representation from gradients
            
            if ~isempty(Grad)
                this.Critic = optimize(this.Critic, Grad.Critic);
                updateTargetRepresentations(this);
            end
        end
        
        function preSettings = preTrain(this)
            preSettings = preTrain@rl.agent.AbstractAgent(this);
            preSettings.DoValidateForExperienceBuffer = ...
                this.ExperienceBuffer.DoValidate;
            this.ExperienceBuffer.DoValidate = false;
        end
        
        function postTrain(this,preSettings)
            postTrain@rl.agent.AbstractAgent(this,preSettings);
            this.ExperienceBuffer.DoValidate = ...
                preSettings.DoValidateForExperienceBuffer;
        end
        
        function actor = getActor(this) %#ok<MANU>
            % getActor: DRQN agent does not have actor. Therefore, it
            % returns empty.
            
            actor = [];
        end
        
        function appendExperience(this,experiences)
            % append experiences to buffer
            append(this.ExperienceBuffer,{experiences});
        end
    end
    
    methods (Static)
        function obj = loadobj(s)
            if isstruct(s)
                if ~isfield(s,'Version')
                    % version 1 does not have Version field
                    % In version 2,
                    %   - Critic changes from rlRepresentation to rlQValueRepresentation
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
                end
                                
                obj = rl.agent.rlDQNAgent(s.Critic,s.AgentOptions_);
                % REVISIT: how to retain the state of optimizer
                obj.ExplorationModel = s.ExplorationModel;
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