function Options = rlPGAgentOptions(varargin)
% rlPGAgentOptions: Creates options for PG Agent.
%
%   OPT = rlPGAgentOptions returns the default options for rlPGAgent. 
%
%   OPT = rlPGAgentOptions('Option1',Value1,'Option2',Value2,...) uses name/value
%   pairs to override the default values for 'Option1','Option2',...
%
%   Supported options include:
%
%   UseBaseline                     Specifies whether to use value function as baseline for learning
%   EntropyLossWeight               Weight for entropy loss to promote policy exploration.
%   SampleTime                      Sample time of the agent
%   DiscountFactor                  Discount factor to apply to future rewards during training
%
%   See also: rlPGAgent, rlDDPGAgentOptions, rlACAgentOptions, rlDQNAgentOptions

% Copyright 2017-2018 The MathWorks Inc.

Options = rl.option.rlPGAgentOptions(varargin{:});

end