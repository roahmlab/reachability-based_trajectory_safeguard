 classdef rlDQNAgentOptions < rl.option.AgentMemoryTarget
    % rlDQNAgentOptions: Creates agent options for DQN
    
    % Copyright 2017-2018 The MathWorks Inc.
    
    % Note: create mixins for different categories    
    properties        
        % Use Double DQN for learning
        UseDoubleDQN
        
        % Options for epsilon greedy exploration
        EpsilonGreedyExploration
        
        % Maximum trajectory length for training for LSTM
        SequenceLength = 1
    end
    
    methods
        function obj = rlDQNAgentOptions(varargin)
            obj = obj@rl.option.AgentMemoryTarget(varargin{:});
            
            parser = obj.Parser;
            % Double DQN
            addParameter(parser,'UseDoubleDQN',true);
            
            % Epsilon Greedy
            addParameter(parser,'EpsilonGreedyExploration',rl.option.EpsilonGreedyExploration),...

            % Maximum trajectory length for training (LSTM)
            addParameter(parser,'SequenceLength',1);

            parse(parser,varargin{:});
            obj.Parser = parser;
            obj.UseDoubleDQN = parser.Results.UseDoubleDQN;
            obj.EpsilonGreedyExploration = parser.Results.EpsilonGreedyExploration;
            obj.SequenceLength = parser.Results.SequenceLength;
            
            parser.KeepUnmatched = false;
            parse(parser,varargin{:})            
        end
        function obj = set.UseDoubleDQN(obj,Value)
            validateattributes(Value,{'logical','numeric'},{'scalar','real','finite'},'','UseDoubleDQN');
            obj.UseDoubleDQN = logical(Value);
        end
        function obj = set.EpsilonGreedyExploration(obj,Value)
            validateattributes(Value,{'rl.option.EpsilonGreedyExploration'},{'scalar'},'','EpsilonGreedyExploration');
            obj.EpsilonGreedyExploration = Value;
        end        
        function obj = set.SequenceLength(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','integer','positive','finite'},'','SequenceLength');
            obj.SequenceLength = Value;
        end
    end
end