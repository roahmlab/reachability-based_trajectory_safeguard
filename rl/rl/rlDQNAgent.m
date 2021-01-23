function Agent = rlDQNAgent(Critic, varargin)
% rlDQNAgent: Creates DQN agent.
%
%   agent = rlDQNAgent(CRITIC) creates a deep Q-network agent with default
%   options and the specified critic representation.
%
%   agent = rlDQNAgent(CRITIC,OPTIONS) creates a deep Q-network agent with
%   the specified options. To create OPTIONS, use rlDQNAgentOptions.
%
%   See also: rlDDPGAgent, rlPGAgent, rlACAgent, rlDQNAgentOptions

% Copyright 2018 The MathWorks, Inc.

narginchk(1,2)
if isempty(varargin)
    Option = rlDQNAgentOptions;
else
    Option = varargin{1};
    validateattributes(Option, {'rl.option.rlDQNAgentOptions'}, {'scalar', 'nonempty'}, mfilename, 'Options', 2);
end
validateattributes(Critic, {'rl.representation.rlQValueRepresentation'}, {'scalar', 'nonempty'}, mfilename, 'Critic', 1);

Agent = rl.agent.rlDQNAgent(Critic, Option);

end