function Options = rlACAgentOptions(varargin)
% rlACAgentOptions: Creates options for AC Agent.
%
%   OPT = rlACAgentOptions returns the default options for rlACAgent. 
%
%   OPT = rlACAgentOptions('Option1',Value1,'Option2',Value2,...) uses name/value
%   pairs to override the default values for 'Option1','Option2',...
%
%   Supported options are:
%
%   NumStepsToLookAhead             Number of steps to look-ahead when computing return
%   EntropyLossWeight               Weight for entropy loss to promote policy exploration.
%   SampleTime                      Sample time of the agent
%   DiscountFactor                  Discount factor to apply to future rewards during training
%
%   See also: rlACAgent, rlDDPGAgentOptions, rlPGAgentOptions, rlDQNAgentOptions

% Copyright 2017-2018 The MathWorks, Inc.

Options = rl.option.rlACAgentOptions(varargin{:});

end