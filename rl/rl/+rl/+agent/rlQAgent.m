classdef rlQAgent < rl.agent.AbstractAgent
    % rlQAgent: Creates a Q-learning agent.
    %
    %   agent = rlQAgent(CRITIC) creates a Q-learning agent with default
    %   options and the specified critic representation.
    %
    %   agent = rlQAgent(CRITIC,OPTIONS) creates a Q-learning agent with
    %   the specified options. To create OPTIONS, use rlQAgentOptions.
    %
    %   See also: rlQAgentOptions, rlSARSAAgent, rlDDPGAgent, rlPGAgent, rlACAgent, rlDQNAgent
    
    % Copyright 2017-2019 The MathWorks Inc.   

    properties (Dependent)
        % Options to configure RL agent
        AgentOptions
    end
    
    properties (Access = private)
        % Private options to configure RL agent
        AgentOptions_ = [];
    end
    
    properties (Access = private)
        % Critic Representation
        Critic
        
        % Epsilon-Greedy Exploration parameters
        ExplorationModel
    end
    
    properties (Constant, Access = private)
        % version indicator for backward compatibility
        Version = 2
    end
    
    methods
        function this = rlQAgent(Critic, Option)
            % Constructor
            
            this = this@rl.agent.AbstractAgent();
            
            % extract observation and action info
            this.ActionInfo = Critic.ActionInfo;
            this.ObservationInfo = Critic.ObservationInfo;
            
            % set agent option
            this.AgentOptions = Option;
            
            % REVISIT: not support agent has state
            if hasState(Critic)
                error(message('rl:agent:errAgentHasStateNotSupport'))
            end
            
            % set representation
            this = setCritic(this,Critic);
            this.HasCritic = true;
  
            % Construct exploration model
            this.ExplorationModel = this.AgentOptions.EpsilonGreedyExploration;
        end
        
        function Action = getActionWithExploration(this,Observation)
            % Return an action with exploration
            
            if rand < this.ExplorationModel.Epsilon
                Action = usample(this.ActionInfo);
            else
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
        function Critic = getCritic(this)
            % getCritic: Return the critic representation, CRITIC, for the
            % specified reinforcement learning agent, AGENT.
            %
            %   CRITIC = getCritic(AGENT)
            
            Critic = this.Critic;
        end
        
        function this = setCritic(this,Critic)
            % setCritic: Set the critic of the reinforcement learning agent
            % using the specified representation, CRITIC, which must be
            % consistent with the observations and actions of the agent.
            %
            %   AGENT = setCritic(AGENT,CRITIC)
            %
                        
            % validate critic is a Q representation
            validateattributes(Critic, {'rl.representation.rlQValueRepresentation'}, {'scalar', 'nonempty'}, '', 'Critic');
            
            % validate action ans observation infos are same
            if ~isCompatible(this.ActionInfo,Critic.ActionInfo)
                error(message('rl:agent:errIncompatibleActionInfo'))
            end
            if ~isCompatible(this.ObservationInfo,Critic.ObservationInfo)
                error(message('rl:agent:errIncompatibleObservationInfo'))
            end
            
            % set loss
            Critic = setLoss(Critic,"rl.loss.rlmse");
            
            % DQN does not support critic with continuous action data spec
            if isa(Critic.ActionInfo, 'rl.util.rlNumericSpec')
                error(message('rl:agent:errDQNContinuousCritic'))
            end
            
            % Set critic network
            this.Critic = Critic;
            
            reset(this);
        end
        
        function set.AgentOptions(this,NewOptions)
            validateattributes(NewOptions,{'rl.option.rlQAgentOptions'},{'scalar'},'','AgentOptions');
            % validate experience buffer has sufficient length
          
            this.AgentOptions_ = NewOptions;
            this.SampleTime = NewOptions.SampleTime;
            this.ExplorationModel = NewOptions.EpsilonGreedyExploration;
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
            % Q agent code gen
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
            q0 = getMaxQValue(this.Critic,observation);
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
            [~,ActionIdx] = getMaxQValue(this.Critic, Observation);
            Action = getElementValue(this.ActionInfo,ActionIdx);
        end
        
        function resetImpl(this)
            % Revert exploration model to original parameters
            this.ExplorationModel = this.AgentOptions.EpsilonGreedyExploration;
        end
        
        function action = learn(this,exp)
            % learn from current experiences, return action with exploration
            % exp = {state,action,reward,nextstate,isdone}
            
            Observation     = exp{1};
            Action          = exp{2};
            Reward          = exp{3};
            NextObservation = exp{4};
            Done            = exp{5};
            
            if Done == 1
                % for final step, just use the immediate reward, since
                % there is no more a next state
                QTargetEstimate = Reward;
            else
                TargetQValue = getMaxQValue(this.Critic, NextObservation);
                QTargetEstimate = Reward + this.AgentOptions.DiscountFactor * TargetQValue;
            end
            
            % update the critic with computed target
            CriticGradient = gradient(this.Critic,'loss-parameters',...
                    [Observation,Action], QTargetEstimate);
            this.Critic = optimize(this.Critic, CriticGradient);
            
            % update exploration model
            this.ExplorationModel = update(this.ExplorationModel);
            
            % compute action from the current policy
            action = getActionWithExploration(this, NextObservation);
        end
        
        function trainingOptions = validateAgentTrainingCompatibilityImpl(~,trainingOptions)
            % Validate Q agent training options compatibility
            
            if ~strcmpi(trainingOptions.Parallelization,'none')
                % currently do not support parallel training
                error(message('rl:general:errParallelTrainNotSupport'));
            end
        end
        
        function HasState = hasStateImpl(this)
            HasState = hasState(this.Critic);
        end
    end
    
    methods (Hidden)
        function actor = getActor(~)
            % getActor: Q agent does not have actor. Therefore, it
            % returns empty.            
            
            actor = [];
        end
    end
    
    methods (Static)
        function obj = loadobj(s)
            if isstruct(s)
                if ~isfield(s,'Version')
                    % version 1 does not have Version field
                    % In version 2,
                    %   - Critic changes from rlRepresentation to rlQValueRepresentation
                    criticModel = s.Critic.getModel;
                    if isa(criticModel,'DAGNetwork') || isa(criticModel,'nnet.cnn.LayerGraph')
                        s.Critic = rlQValueRepresentation(criticModel,...
                            s.ObservationInfo, s.ActionInfo, ...
                            'Observation', s.Critic.ObservationNames, ...
                            'Action', s.Critic.ActionNames, ...
                            s.Critic.Options);
                    else
                        s.Critic = rlQValueRepresentation(criticModel,...
                            s.ObservationInfo, s.ActionInfo, ...
                            s.Critic.Options);
                    end
                end
                obj = rl.agent.rlQAgent(s.Critic,s.AgentOptions_);
                % REVISIT: how to retain the state of optimizer
                obj.ExplorationModel = s.ExplorationModel;
            else
                obj = s;
            end
        end
    end
end
