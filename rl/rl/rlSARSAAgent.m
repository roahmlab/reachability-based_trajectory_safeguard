function Agent = rlSARSAAgent(Critic, varargin)
% rlQAgent: Creates a SARSA agent.
%
%   agent = rlSARSAAgent(CRITIC) creates a SARSA agent with default
%   options and the specified critic representation.
%
%   agent = rlSARSAAgent(CRITIC,OPTIONS) creates a SARSA agent with
%   the specified options. To create OPTIONS, use rlSARSAAgentOptions.
%
%   See also: rlSARSAAgentOptions, rlQAgent, rlDDPGAgent, rlPGAgent, rlACAgent, rlDQNAgent

% Copyright 2017-2019 The MathWorks Inc.

narginchk(1,2)
if isempty(varargin)
    Option = rlSARSAAgentOptions;
else
    Option = varargin{1};
    validateattributes(Option, {'rl.option.rlSARSAAgentOptions'}, {'scalar', 'nonempty'}, mfilename, 'Options', 2);
end
validateattributes(Critic, {'rl.representation.rlQValueRepresentation'}, {'scalar', 'nonempty'}, mfilename, 'Critic', 1);

Agent = rl.agent.rlSARSAAgent(Critic, Option);

end