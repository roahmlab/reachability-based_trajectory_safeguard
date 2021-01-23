classdef EpsilonGreedyExploration
    % EpsilonGreedyExploration: Creates options for Epsilon-greedy
    % exploration for DQN agents. You can also create EpsilonGreedyExploration
    % objects for your own custom agents.
    %
    % Supported Options are:
    %
    %     Epsilon: Probability threshold for the agent to either randomly
    %     select a valid action or select the action that maximizes the
    %     state-action value function. A larger value of Epsilon implies
    %     that the probability is larger for the agent to select actions
    %     randomly; in other words, the exploration rate is higher.
    %
    %     EpsilonMin: Threshold that decides whether Epsilon is updated
    %
    %     EpsilonDecay: Decay rate of Epsilon when Epsilon is updated
    %
    % Epsilon is updated according to the following formula when it is
    % greater than EpsilonMin:
    %
    %    Epsilon = Epsilon*(1-EpsilonDecay)
    %
    % Examples:
    %
    % obj = rl.option.EpsilonGreedyExploration;
    % obj = rl.option.EpsilonGreedyExploration('Epsilon',0.8);
    % obj = rl.option.EpsilonGreedyExploration('Epsilon',0.8,...
    %                                          'EpsilonMin',0.2);
    % obj = rl.option.EpsilonGreedyExploration('Epsilon',0.8,...
    %                                          'EpsilonMin',0.02,...
    %                                          'EpsilonDecay',1e-3);
    % obj = update(obj); % updates Epsilon
    %
    % An EpsilonGreedyExploration object is automatically created within 
    % an rlDQNAgentOptions object. In this case, specify 
    % EpsilonGreedyExploration options using dot notation.
    %
    % opt = rlDQNAgentOptions;
    % opt.EpsilonGreedyExploration.Epsilon = 0.8;

    % Copyright 2019 The MathWorks Inc.
    
    properties
        % Factor by which the exploration rate Epsilon is
        % adjusted each time the agent replays training from its replay
        % buffer (0: no decay, 1: immediate decay)
        EpsilonDecay
    end
    
    properties (Dependent)
        % Rate of exploration vs using the model for
        % predicting the next state during training. Set it to 1 if you
        % want to take random actions at the beginning of training.
        Epsilon
        
        % Final rate of exploration. Set it as low
        % as possible to avoid taking random actions as the model improves
        % or trains on data generated during the training process
        EpsilonMin
    end
    
    properties (Access = private)
        % Rate of exploration vs using the model for
        % predicting the next state during training. Set it to 1 if you
        % want to take random actions at the beginning of training.
        Epsilon_
        
        % Final rate of exploration. Set it as low
        % as possible to avoid taking random actions as the model improves
        % or trains on data generated during the training process
        EpsilonMin_
    end
    
    methods
        function obj = EpsilonGreedyExploration(varargin)
            parser = inputParser();            
            addParameter(parser,'Epsilon',1)
            addParameter(parser,'EpsilonMin',0.01)
            addParameter(parser,'EpsilonDecay',0.005)
            parse(parser,varargin{:}); 
            
            obj.Epsilon = parser.Results.Epsilon;
            obj.EpsilonMin = parser.Results.EpsilonMin;
            obj.EpsilonDecay = parser.Results.EpsilonDecay;
            
            % Following condition guards against Epsilon < EpsilonMin after
            % object construction
            if obj.Epsilon_ < obj.EpsilonMin_
                warning(message('rl:agent:warnEpsilonLtEpsilonMin'));
            end
        end
        
        function obj = set.Epsilon(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','positive','>',0,'<=',1},'','Epsilon');
            % isempty condition required to avoid warnings while loading
            % objects saved with Epsilon < EpsilonMin
            if ~isempty(obj.EpsilonMin_) && ~isempty(obj.Epsilon_)
                if Value < obj.EpsilonMin_
                    warning(message('rl:agent:warnEpsilonLtEpsilonMin'));
                end
            end
            obj.Epsilon_ = Value;
        end
        
        function obj = set.EpsilonMin(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','nonnegative','>',0,'<=',1},'','EpsilonMin');
            % isempty condition required to avoid warnings while loading
            % objects saved with Epsilon < EpsilonMin
            if ~isempty(obj.Epsilon_) && ~isempty(obj.EpsilonMin_)
                if Value > obj.Epsilon_
                    warning(message('rl:agent:warnEpsilonLtEpsilonMin'));
                end
            end
            obj.EpsilonMin_ = Value;
        end
        
        function obj = set.EpsilonDecay(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','positive','>',0,'<=',1},'','EpsilonDecay')
            obj.EpsilonDecay = Value;
        end        
        
        function obj = update(obj)
            obj.Epsilon = max(obj.Epsilon * (1-obj.EpsilonDecay), obj.EpsilonMin);
        end
        
        function Options = get.Epsilon(this)
            Options = this.Epsilon_;
        end
        
        function Options = get.EpsilonMin(this)
            Options = this.EpsilonMin_;
        end
    end
end