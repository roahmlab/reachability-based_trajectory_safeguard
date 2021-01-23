classdef CustomAgent < rl.agent.AbstractAgent
    % CUSTOMAGENT
    %   Subclass from this class when creating custom agents. See the
    %   LQRCustomAgent example for more details.
    
    % Revised: 2-13-2019 
    
    % Copyright 2019 The MathWorks, Inc.
    
    %% methods required for implementation
    methods (Abstract,Access = protected)
        % learnImpl
        %   Define how the agent learns from an experience (exp). exp is a
        %   cell array with the following format:
        %       exp = {observation,action,reward,nextObservation,isDone}
        action = learnImpl(this,exp)
        % getActionWithExplorationImpl
        %   Define the agent exploration. For agents like DDPG, action =
        %   getAction(observation) is called and noise is added to action.
        action = getActionWithExplorationImpl(this,observation)
    end
    % methods to be implemented from AbstractPolicy:
    %   getActionImpl
    %       Define how an action is generated from an observation.
    %   action = getActionImpl(this,observation)
    
    %% default implementations
    methods (Access = protected)
        %   resetImpl(this)
        %       Define how the agent is reset before training.
        function resetImpl(this) %#ok<MANU>
            % Overload for custom behavior
        end
        
        % Validate agent training options compatibility
        function trainingOptions = validateAgentTrainingCompatibilityImpl(~,trainingOptions)
            if ~strcmpi(trainingOptions.Parallelization,'none')
                % currently do not support parallel training
                error(message('rl:general:errParallelTrainNotSupport'));
            end
        end
        
        % set/get tunable parameters
        function setLearnableParametersImpl(this,p) %#ok<INUSD>
            % null op
        end
        function p = getLearnableParametersImpl(this) %#ok<MANU>
            % return []
            p = [];
        end
        
        function HasState = hasStateImpl(this) %#ok<MANU>
            % HASSTATE(AGENT) indicates whether the agent is stateful
            % (use recurrent neural networks) or not. Agent subclass must
            % overload this method for stateful training.
            
            HasState = false;
        end
    end
    %% implementations
    methods (Access = protected)
        function action = learn(obj,exp)
            action = learnImpl(obj,exp);
        end
    end
    
    methods
        function action = getActionWithExploration(this,observation)
            action = getActionWithExplorationImpl(this,observation);
        end
    end
end