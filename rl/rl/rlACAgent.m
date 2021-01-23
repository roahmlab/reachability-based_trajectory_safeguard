function Agent = rlACAgent(Actor, Critic, varargin)
% rlACAgent: Creates AC agent.
%
%   agent = rlACAgent(ACTOR,CRITIC) creates an actor-critic agent with 
%   default options using specified actor and critic networks.
%
%   agent = rlACAgent(ACTOR,CRITIC,OPTIONS) creates an actor-critic agent
%   with the specified options. To create OPTIONS, use rlACAgentOptions. 
%
%   See also: rlDQNAgent, rlDDPGAgent, rlPGAgent, rlACAgentOptions

% Copyright 2018 The MathWorks, Inc.

narginchk(2,3)
if isempty(varargin)
    Options = rlACAgentOptions;
else
    Options = varargin{1};
    validateattributes(Options, {'rl.option.rlACAgentOptions'}, {'scalar', 'nonempty'}, mfilename, 'Options', 3);
end

validateattributes(Actor, {'rl.representation.rlStochasticActorRepresentation','rl.util.rlAbstractRepresentation'}, {'scalar', 'nonempty'}, mfilename, 'Actor', 1);
validateattributes(Critic, {'rl.representation.rlValueRepresentation','rl.util.rlAbstractRepresentation'}, {'scalar', 'nonempty'}, mfilename, 'Critic', 2);

Agent = rl.agent.rlACAgent(Actor, Critic, Options);

end