function agent = rlSACAgentNoTanh(Actor, Critic, varargin)
% rlSACAgent: Creates a SAC agent.
%
%   agent = rlSACAgentNoTanh(ACTOR,CRITIC) creates a soft actor critic agent
%   with default options using the specified actor
%   and critic representations.
%
%   agent = rlSACAgentNoTanh(ACTOR,CRITIC,OPTIONS) creates a soft
%   actor critic agent with the specified options. To
%   create OPTIONS, use rlSACAgentOptions.
%
%   See also: rlSACAgentOptions, EntropyWeightOptions, rlDQNAgent, rlPGAgent, rlACAgent, 
 
% Copyright 2019 The MathWorks, Inc.

narginchk(2,4)
validateattributes(Actor,{'rl.representation.rlStochasticActorRepresentation'},{'scalar','nonempty'},mfilename,'Actor',1);
validateattributes(Critic,{'rl.representation.rlQValueRepresentation'},{'scalar','nonempty'},mfilename,'Critic1',2);

agent = rl.agent.rlSACAgentNoTanh(Actor, Critic, varargin{:});
end