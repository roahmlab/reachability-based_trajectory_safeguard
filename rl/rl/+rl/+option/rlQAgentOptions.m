classdef rlQAgentOptions < rl.option.AgentGeneric
    % rlQAgentOptions: Creates options for Q Learning Agent.
    %
    %   OPT = rlQAgentOptions returns the default options for rlQAgent.
    %
    %   OPT = rlQAgentOptions('Option1',Value1,'Option2',Value2,...) uses name/value
    %   pairs to override the default values for 'Option1','Option2',...
    %
    %   Supported options are:
    %
    %   EpsilonGreedyExploration            Parameters for Epsilon Greedy exploration
    %       Epsilon                         Probability threshold for agent to either randomly
    %                                       select a valid action or select the action that
    %                                       maximizes the state-action value function
    %       EpsilonMin                      Minimum value of Epsilon
    %       EpsilonDecay                    Decay rate of Epsilon when Epsilon is updated
    %   SampleTime                          Sample time of the agent
    %   DiscountFactor                      Discount factor to apply to future rewards during training
    %
    %   See also: rlQAgent, rlDDPGAgentOptions, rlPGAgentOptions, rlACAgentOptions
    
    % Copyright 2017-2018 The MathWorks Inc.
    
    properties               
        % Options for epsilon greedy exploration
        EpsilonGreedyExploration
    end
    
    methods
        function obj = rlQAgentOptions(varargin)
            obj = obj@rl.option.AgentGeneric(varargin{:});
            
            parser = obj.Parser;
           
            % Epsilon Greedy
            addParameter(parser,'EpsilonGreedyExploration',rl.option.EpsilonGreedyExploration),...
            
            parse(parser,varargin{:});
            obj.Parser = parser;
            obj.EpsilonGreedyExploration = parser.Results.EpsilonGreedyExploration;
            
            parser.KeepUnmatched = false;
            parse(parser,varargin{:})            
        end

        function obj = set.EpsilonGreedyExploration(obj,Value)
            validateattributes(Value,{'rl.option.EpsilonGreedyExploration'},{'scalar'},'','EpsilonGreedyExploration');
            obj.EpsilonGreedyExploration = Value;
        end        
    end
end