function Agent = rlQAgent(Critic, varargin)
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

narginchk(1,2)
if isempty(varargin)
    Option = rlQAgentOptions;
else
    Option = varargin{1};
    validateattributes(Option, {'rl.option.rlQAgentOptions'}, {'scalar', 'nonempty'}, mfilename, 'Options', 2);
end
validateattributes(Critic, {'rl.representation.rlQValueRepresentation'}, {'scalar', 'nonempty'}, mfilename, 'Critic', 1);

Agent = rl.agent.rlQAgent(Critic, Option);

end