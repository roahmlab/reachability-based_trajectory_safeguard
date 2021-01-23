function agent = rlTD3Agent(Actor, Critic, varargin)
% RLTD3AGENT: Creates TD3 agent.
%
%   AGENT = RLTD3AGENT(ACTOR,CRITIC) creates a delayed deep deterministic 
%   (delayed DDPG) agent with specified actor and critic represenatations 
%   with default options. It is a DDPG agent with target policy smoothing
%   and delayed policy and targets update.
%       ACTOR: rlDeterministicActorRepresentation
%       CRITIC1: rlQValueRepresentation
%
%   AGENT = RLTD3AGENT(ACTOR,CRITICS) creates a twin delayed deep
%   deterministic policy gradient agent with default options using the
%   specified ACTOR and 1x2 vector CRITICS representations. Each critic in
%   CRITICS vector must have different estimation for the same observation 
%   and action pair. For example, they can have same structure but
%   different initial parameters or have different structure.
%
%   For all of the previous syntaxes, you can specify nondefault options
%   using an rlTD3AgentOptions object, OPTIONS.
% 
%       AGENT = RLTD3AGENT(...,OPTIONS)
%
%   See also: rlDQNAgent, rlPPOAgent, rlDDPGAgent, rlTD3AgentOptions
 
% Copyright 2019 The MathWorks, Inc.

narginchk(2,3)
if isempty(varargin)
    Option = rlTD3AgentOptions;
else
    Option = varargin{1};
    validateattributes(Option, {'rl.option.rlTD3AgentOptions'}, {'scalar', 'nonempty'}, mfilename, 'Options', 3);
end

validateattributes(Actor, {'rl.representation.rlDeterministicActorRepresentation'}, {'scalar', 'nonempty'}, mfilename, 'Actor', 1);
validateattributes(Critic, {'rl.representation.rlQValueRepresentation'}, {'vector', 'nonempty'}, mfilename, 'Critic', 2);
agent = rl.agent.rlTD3Agent(Actor, Critic, Option);

