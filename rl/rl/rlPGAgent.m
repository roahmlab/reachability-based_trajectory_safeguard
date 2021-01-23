function Agent = rlPGAgent(Actor, varargin)
% rlPGAgent: Creates PG agent.
%
%   agent = rlPGAgent(ACTOR) a policy gradient agent with default options
%   using the specified actor representation and no critic representation.
%
%   agent = rlPGAgent(ACTOR,OPTIONS) creates a policy gradient agent only
%   agent with actor and options. To create OPTIONS, use rlPGAgentOptions.
%
%   agent = rlPGAgent(ACTOR,CRITIC) creates a policy gradient network with
%   default options using the specified actor and critic representations.
%
%   agent = rlPGAgent(ACTOR,CRITIC,OPTIONS) creates an actor-critic agent
%   with the specified actor, critic, and options. To create OPTIONS, use
%   rlPGAgentOptions.
%
%   See also: rlDQNAgent, rlDDPGAgent, rlPGAgentOptions

% Copyright 2017-2018 The MathWorks Inc.

narginchk(1,3)
validateattributes(Actor, {'rl.representation.rlStochasticActorRepresentation','rl.util.rlAbstractRepresentation'}, {'scalar', 'nonempty'}, mfilename, 'Actor', 1);

switch nargin
    case 1
        % not use base line and use default options
        Critic = [];
        Options = rlPGAgentOptions;
        Options.UseBaseline = false;
    case 2
        validateattributes(varargin{1}, {'rl.representation.rlValueRepresentation','rl.option.rlPGAgentOptions'}, {'scalar', 'nonempty'}, mfilename, '', 2);
        if isa(varargin{1},'rl.option.rlPGAgentOptions')
            % specify options and not use baseline
            Critic = [];
            Options = varargin{1};
            if Options.UseBaseline
                error(message('rl:agent:errPGBaselineTrueNoCritic'))
            end
        elseif isa(varargin{1},'rl.representation.rlValueRepresentation')
            % use baseline with default options
            Critic = varargin{1};
            Options = rlPGAgentOptions;
        end
    case 3
        % specify both baseline and options
        Critic = varargin{1};
        Options = varargin{2};
        validateattributes(Critic, {'rl.representation.rlValueRepresentation'}, {'scalar', 'nonempty'}, mfilename, 'Critic', 2);
        validateattributes(Options, {'rl.option.rlPGAgentOptions'}, {'scalar', 'nonempty'}, mfilename, 'Options', 3);
        if ~Options.UseBaseline
            warning(message('rl:agent:warnPGBaselineFalseHaveCritic'))
        end
end

Agent = rl.agent.rlPGAgent(Actor, Critic, Options);

end