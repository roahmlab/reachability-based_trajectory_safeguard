function env = BasicGridWorld 
% BASICGRIDWORLD: Creates basic gridworld RL environment.

% Copyright 2018 The MathWorks, Inc.

% construct basic gridworld environment
GW = createGridWorld(5,5);

% set current, terminal states and obstacle states
GW.TerminalStates = "[5,5]";
GW.ObstacleStates = ["[3,3]";"[3,4]";"[3,5]";"[4,3]"];

% update gridworld
updateStateTranstionForObstacles(GW)

% set jump rule
GW.T(state2idx(GW,"[2,4]"),:,:) = 0;
GW.T(state2idx(GW,"[2,4]"),state2idx(GW,"[4,4]"),:) = 1;

% set rewards
nS = numel(GW.States);
nA = numel(GW.Actions);
GW.R = -1*ones(nS,nS,nA);
GW.R(state2idx(GW,"[2,4]"),state2idx(GW,"[4,4]"),:) = 5;
GW.R(:,state2idx(GW,GW.TerminalStates),:) = 10;

env = rlMDPEnv(GW);


    
    
            

    
