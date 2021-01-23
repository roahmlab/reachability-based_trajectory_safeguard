function Options = rlQAgentOptions(varargin)
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

Options = rl.option.rlQAgentOptions(varargin{:});

end