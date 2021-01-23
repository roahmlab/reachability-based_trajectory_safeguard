classdef AbstractOnPolicyPGAgent < rl.agent.AbstractAgent
    %ABSTRACTONPOLICYPGAGENT Parent class of on-policy policy gradient
    % agent (PG, PPO, AC).
    
    % Copyright 2019 The MathWorks Inc.
    
    properties (Dependent)
        % Options to configure RL agent
        AgentOptions
    end
    
    properties (Access = protected)
        % Private options to configure RL agent
        AgentOptions_ = [];
        
        % Critic function approximator
        Critic
        
        % Actor function approximator
        Actor
    end
    
    properties (SetAccess = protected, Hidden)
        ExperienceBuffer
    end
    
    methods
        function this = AbstractOnPolicyPGAgent(Actor, Critic, Options)
            % Constructor
            
            this = this@rl.agent.AbstractAgent();
            
            % extract observation and action info
            this.ObservationInfo = Actor.ObservationInfo;
            this.ActionInfo = Actor.ActionInfo;
            
            % set options
            this.AgentOptions = Options;
            
            % set representations
            setActor(this, Actor);
            setCritic(this, Critic);
        end
        
        function Action = getActionWithExploration(this,Observation)
            % GETACTIONWITHEXPLORATION return actions sampled from the
            % probability distribution of the actor policy given  
            % observations.
            %
            %   ACTION = GETACTIONWITHEXPLORATION(AGENT, OBSERVATION) 
            %   return ACTION sampled from the current probability
            %   distribution of the actor policy of AGENT given OBSERVATION.
            
            Action = getActionImpl(this,Observation);
        end
    end
    
    %======================================================================
    % Implementation of abstract methods
    %======================================================================
    methods (Access = protected)
        % policy method
        function Action = getActionImpl(this,Observation)
            % Given the current observation of the system, return an
            % action with highest probability from the actor policy
            
            if hasState(this)
                % reset hidden state of actor at the beginning of each
                % train/sim episode
                if this.EpisodeInfo.StepsTaken == 0
                    this.Actor = resetState(this.Actor);
                end
            end
            [Action, State] = getAction(this.Actor, Observation);
            this.Actor = setState(this.Actor, State);
            
            % return numeric if action is single channel for backward
            % compatibility
            if numel(Action) == 1
                Action = Action{1};
            end
        end
        
        % critic method
        function q0 = evaluateQ0(this,exp)
            % Get the estimated long-term return (Q0) from the agent's
            % critic
            
            if this.HasCritic
                observation = exp{1};
                if hasState(this)
                    % reset hidden state before estimate Q at step 0
                    this.Critic = resetState(this.Critic);
                end
                q0 = getValue(this.Critic,observation);
                if isa(q0,'dlarray')
                    q0 = extractdata(q0);
                end
            else
                q0 = 0;
            end
        end
        
        % stateful check
        function HasState = hasStateImpl(this)
            % all on policy agent has an actor
            HasState = hasState(this.Actor);
        end
        
        % codegen
        function [rep,argStruct] = generateProcessFunctions_(this,argStruct)
            rep = this.Actor;
            if rl.util.isaSpecType(this.ActionInfo,'discrete')
                argStruct = rl.codegen.generateDiscreteStochasticPolicyFcn(argStruct,this.ActionInfo);
            elseif rl.util.isaSpecType(this.ActionInfo,'continuous')
                argStruct = rl.codegen.generateContinuousStochasticPolicyFcn(argStruct,this.ActionInfo);
            else
                error(message('rl:agent:CodegenNotSupportMixedActionType'))
            end
        end
        
        % set/get tunable parameters
        function setLearnableParametersImpl(this,p)
            this.Actor  = setLearnableParameters(this.Actor ,p.Actor );
            if this.HasCritic
                this.Critic = setLearnableParameters(this.Critic,p.Critic);
            end
        end
        function p = getLearnableParametersImpl(this)
            p.Actor  = getLearnableParameters(this.Actor );
            if this.HasCritic
                p.Critic = getLearnableParameters(this.Critic);
            end
        end
        
        % copy
        function that = copyElement(this)
            that = copyElement@rl.agent.AbstractAgent(this);
            that.ExperienceBuffer = copy(this.ExperienceBuffer);
        end
    end
    
    %======================================================================
    % Get/set
    %======================================================================
    methods
        function this = setCritic(this, Critic)
            % setCritic: Set the critic of the reinforcement learning agent
            % using the specified representation, CRITIC, which must be
            % consistent with the observations and actions of the agent.
            %
            %   AGENT = setCritic(AGENT,CRITIC)
            
            if isempty(Critic)
                this.Critic = [];
            else
                % validate critic is a value representation
                validateattributes(Critic, {'rl.representation.rlValueRepresentation'}, {'scalar', 'nonempty'}, '', 'Critic');
                
                % validate observation infos are same
                if ~isCompatible(this.ObservationInfo,Critic.ObservationInfo)
                    error(message('rl:agent:errIncompatibleObservationInfo'))
                end
                
                % check if actor and critic created from same data specs and if
                % both stateless or both have state
                rl.agent.AbstractAgent.validateActorCriticInfo(Critic,this.Actor)
                
                % set mse loss
                Critic = setLoss(Critic,"rl.loss.rlmse");
                
                % Set critic representation
                this.Critic = Critic;
            end
            
            % Reset experience buffer (on-policy needs fresh data)
            reset(this);
            
            this.HasCritic = ~isempty(this.Critic);
        end
        
        function Critic = getCritic(this)
            % getCritic: Return the critic representation, CRITIC, for the
            % specified reinforcement learning agent, AGENT.
            %
            %   CRITIC = getCritic(AGENT)
            
            if this.HasCritic
                Critic = this.Critic;
            else
                Critic = [];
            end
        end
        
        function this = setActor(this, Actor)
            % setActor: Set the actor of the reinforcement learning agent
            % using the specified representation, ACTOR, which must be
            % consistent with the observations and actions of the agent.
            %
            %   AGENT = setActor(AGENT,ACTOR)
            
            % validate actor is a stochastic actor representation
            validateattributes(Actor, {'rl.representation.rlStochasticActorRepresentation'}, {'scalar', 'nonempty'}, '', 'Actor');
            
            % validate action and observation infos are same
            if ~isCompatible(this.ActionInfo, Actor.ActionInfo)
                error(message('rl:agent:errIncompatibleActionInfo'))
            end
            if ~isCompatible(this.ObservationInfo, Actor.ObservationInfo)
                error(message('rl:agent:errIncompatibleObservationInfo'))
            end
            
            % check if actor and critic created from same data specs and if
            % both stateless or both have state
            rl.agent.AbstractAgent.validateActorCriticInfo(this.Critic,Actor)
            
            % set actor loss
            Actor = setActorLoss(this, Actor);
            
            % Set actor representation
            this.Actor = Actor;
            
            % Reset experience buffer (on-policy needs fresh data)
            reset(this);
        end
        
        function Actor = getActor(this)
            % getActor: Return the actor representation, ACTOR, for the
            % specified reinforcement learning agent.
            %
            %   ACTOR = getActor(AGENT)
            
            Actor = this.Actor;
        end
        
        function set.AgentOptions(this,NewOptions)
            validateOptionImpl(this,NewOptions);
            this.AgentOptions_ = NewOptions;
            this.SampleTime = NewOptions.SampleTime;
        end
        function Options = get.AgentOptions(this)
            Options = this.AgentOptions_;
        end
    end
    
    methods (Hidden)
        function appendExperience(this,experiences)
            % append experiences to buffer
            append(this.ExperienceBuffer,{experiences});
        end
    end
    
    methods (Access = protected)
        validateOptionImpl(this, NewOptions)
        Actor = setActorLoss(this, Actor)
    end
    
    methods (Static)
        function [Actor, Critic] = redirectV1Rep(s)
            if ~isfield(s,'Version')
                % version 1 does not have Version field
                % In version 2,
                %   - Critic changes from rlRepresentation to rlValueRepresentation
                %   - Actor changes from rlRepresentation to rlStochasticActorRepresentation
                %   - ExperienceBuffer requires obs and act dims inputs
                %   but will always reconstruct since agent is on-policy
                
                ActorModel = s.Actor.getModel;
                if isa(ActorModel,'DAGNetwork') || isa(ActorModel,'nnet.cnn.LayerGraph')
                    Actor = rlStochasticActorRepresentation(ActorModel,...
                        s.ObservationInfo, s.ActionInfo, ...
                        'Observation', s.Actor.ObservationNames, ...
                        s.Actor.Options);
                else
                    Actor = rlStochasticActorRepresentation(ActorModel,...
                        s.ObservationInfo, s.ActionInfo, ...
                        s.Actor.Options);
                end
                
                if ~isempty(s.Critic)
                    CriticModel = s.Critic.getModel;
                    if isa(CriticModel,'DAGNetwork') || isa(CriticModel,'nnet.cnn.LayerGraph')
                        Critic = rlValueRepresentation(CriticModel,...
                            s.ObservationInfo, ...
                            'Observation', s.Critic.ObservationNames, ...
                            s.Critic.Options);
                    else
                        Critic = rlValueRepresentation(CriticModel,...
                            s.ObservationInfo, ...
                            s.Critic.Options);
                    end
                else
                    Critic = s.Critic;
                end
            end
        end
    end
end
