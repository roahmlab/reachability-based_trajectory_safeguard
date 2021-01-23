function env = rlMDPEnv(MDP)
% rlMDPEnv: Create a MATLAB based reinforcement learning environment for a 
% MDP(Markov Decision Process) by supplying the MDP model.
%
%   ENV = rlMDPEnv(MDP) creates a reinforcement learning environment with
%   the specified MDP model. See createGridWorld and createMDP on how to
%   create MDP models.
%
%   Example:
%       % Create a MDP: 5x5 Grid World
%       MDP = createGridWorld(5,5);
%       % Set terminal state the "[5,5]" grid
%       MDP.TerminalStates = "[5,5]";
%       % Set the reward to the terminal state to 10 
%       % and -1 for each step taken
%       nS = numel(MDP.States);
%       nA = numel(MDP.Actions);
%       MDP.R = -1*ones(nS,nS,nA);
%       MDP.R(:,state2idx(MDP,MDP.TerminalStates),:) = 10;
%       % Create the MDP reinforcement learning environment
%       env = rlMDPEnv(MDP);
%       % show visualization
%       plot(env)
%       % Take a step in the environment
%       step(env,1); % Move south
%
%   See also: createGridWorld, createMDP

%   Copyright 2017-2019 The MathWorks, Inc.

env = rl.env.rlMDPEnv(MDP);

end 