function GW = createGridWorld(m,n,varargin)
%createGridWorld Creates a grid world markov decision process
%
%   GW = createGridWorld(m,n) creates a grid world of size m by n with
%   default actions of ["N";"S";"E";"W"]
%
%   GW = createGridWorld(m,n,Moves) creates a grid world of size m by n. 
%   For Moves = "Standard" the actions are ["N";"S";"E";"W"].
%   For Moves = "Kings" the actions are ["N";"S";"E";"W";"NE";"NW";"SE";"SW"]
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
%
%   See also: rlMDPEnv, createMDP

% Copyright 2018 The MathWorks, Inc.

narginchk(2,3)

GW = rl.env.GridWorld(m,n,varargin{:});