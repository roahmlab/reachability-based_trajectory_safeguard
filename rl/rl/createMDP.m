function GW = createMDP(m,n)
%createMDP Creates a markov decision process model
%
%   MDP = createMDP(3,2)
%     creates an MDP with 3 states and 2 actions with default
%     state and action names
%
%   MDP = createMDP(["state1";"state2";"state3"],["action1","action2"])
%     creates an MDP with 3 states and 2 actions with the
%     specified state and action names
%
%   Example:
%       % Create MDP with 8 states and 2 actions
%       MDP = createMDP(8,["up";"down"]);
% 
%       % State 1 Transition and Reward
%       MDP.T(1,2,1) = 1;  MDP.R(1,2,1) = 3;
%       MDP.T(1,3,2) = 1;  MDP.R(1,3,2) = 1;
%       % State 2 Transition and Reward
%       MDP.T(2,4,1) = 1;  MDP.R(2,4,1) = 2;
%       MDP.T(2,5,2) = 1;  MDP.R(2,5,2) = 1;
%       % State 3 Transition and Reward
%       MDP.T(3,5,1) = 1;  MDP.R(3,5,1) = 2;
%       MDP.T(3,6,2) = 1;  MDP.R(3,6,2) = 4;
%       % State 4 Transition and Reward
%       MDP.T(4,7,1) = 1;  MDP.R(4,7,1) = 3;
%       MDP.T(4,8,2) = 1;  MDP.R(4,8,2) = 2;
%       % State 5 Transition and Reward
%       MDP.T(5,7,1) = 1;  MDP.R(5,7,1) = 1;
%       MDP.T(5,8,2) = 1;  MDP.R(5,8,2) = 9;
%       % State 6 Transition and Reward
%       MDP.T(6,7,1) = 1;  MDP.R(6,7,1) = 5;
%       MDP.T(6,8,2) = 1;  MDP.R(6,8,2) = 1;
%       % State 7 Transition and Reward
%       MDP.T(7,7,1) = 1;  MDP.R(7,7,1) = 0;
%       MDP.T(7,7,2) = 1;  MDP.R(7,7,2) = 0;
%       % State 8 Transition and Reward
%       MDP.T(8,8,1) = 1;  MDP.R(8,8,1) = 0;
%       MDP.T(8,8,2) = 1;  MDP.R(8,8,2) = 0;
% 
%       MDP.TerminalStates = ["s7";"s8"];
%
%   See also: rlMDPEnv, createGridWorld

% Copyright 2018 The MathWorks, Inc.

narginchk(2,3)

GW = rl.env.GenericMDP(m,n);