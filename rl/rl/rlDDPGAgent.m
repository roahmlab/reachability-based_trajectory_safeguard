function agent = rlDDPGAgent(Actor, Critic, varargin)
% rlDDPGAgent: Creates DDPG agent.
%
%   agent = rlDDPGAgent(ACTOR,CRITIC) creates a deep deterministic
%   policy gradient agent with default options using the specified actor
%   and critic representations.
%
%   agent = rlDDPGAgent(ACTOR,CRITIC,OPTIONS) creates a deep
%   deterministic policy gradient agent with the specified options. To
%   create OPTIONS, use rlDDPGAgentOptions.
%
%   See also: rlDQNAgent, rlPGAgent, rlACAgent, rlDDPGAgentOptions
 
% Copyright 2018 The MathWorks, Inc.

narginchk(2,3)
if isempty(varargin)
    Option = rlDDPGAgentOptions;
else
    Option = varargin{1};
    validateattributes(Option, {'rl.option.rlDDPGAgentOptions'}, {'scalar', 'nonempty'}, mfilename, 'Options', 3);
end
validateattributes(Actor, {'rl.representation.rlDeterministicActorRepresentation'}, {'scalar', 'nonempty'}, mfilename, 'Actor', 1);
validateattributes(Critic, {'rl.representation.rlQValueRepresentation'}, {'scalar', 'nonempty'}, mfilename, 'Critic', 2);

agent = rl.agent.rlDDPGAgent(Actor, Critic, Option);

