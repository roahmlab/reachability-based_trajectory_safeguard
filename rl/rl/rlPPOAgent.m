function Agent = rlPPOAgent(Actor, Critic, varargin)
% RLPPOAGENT: Creates Proximal Policy Optimization (PPO) agent.
%
%   AGENT = RLPPOAGENT(ACTOR,CRITIC) creates a PPO agent with 
%   default options using specified actor and critic networks.
%
%   AGENT = RLPPOAGENT(ACTOR,CRITIC,OPTIONS) creates a PPO agent
%   with the specified options. To create OPTIONS, use rlPPOAgentOptions. 
%
%   See also: rlDQNAgent, rlDDPGAgent, rlPGAgent, rlACAgent, rlPPOAgentOptions

% Copyright 2019 The MathWorks, Inc.

narginchk(2,3)
if isempty(varargin)
    Options = rlPPOAgentOptions;
else
    Options = varargin{1};
    validateattributes(Options, {'rl.option.rlPPOAgentOptions'}, {'scalar', 'nonempty'}, mfilename, 'Options', 3);
end

validateattributes(Actor, {'rl.representation.rlStochasticActorRepresentation','rl.util.rlAbstractRepresentation'}, {'scalar', 'nonempty'}, mfilename, 'Actor', 1);
validateattributes(Critic, {'rl.representation.rlValueRepresentation','rl.util.rlAbstractRepresentation'}, {'scalar', 'nonempty'}, mfilename, 'Critic', 2);

Agent = rl.agent.rlPPOAgent(Actor, Critic, Options);

end