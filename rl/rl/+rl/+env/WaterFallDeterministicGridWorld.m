function env = WaterFallDeterministicGridWorld 
% WaterFallDeterministicGridWorld: Creates Deterministic Water Fall RL environment.

% Copyright 2018 The MathWorks, Inc.

%% Set up basic GridWorld environment
nRow = 8;    % number of rows
nCol = 7;   % number of columns
envWaterDeterministic = createGridWorld(nRow,nCol); 
startState = "[5,1]";
envWaterDeterministic.CurrentState = startState;
envWaterDeterministic.TerminalStates = "[4,5]";
ActionsList = envWaterDeterministic.Actions;        % get list of possible action
StatesList = envWaterDeterministic.States;          % get list of states
stateSize = numel(StatesList);    % number of states
actionSize = numel(ActionsList);   % number of possible actions 

%% Set Rewards Transitions
envWaterDeterministic.R = -1*ones(stateSize,stateSize,actionSize); 
envWaterDeterministic.R(:,state2idx(envWaterDeterministic,"[4,5]"),:) = 10;

%% Set up waterfall rule
% Create waterfall map with water intesity as specified in example 6.5
% From default gridworld, find transitions in waterfall zone and change
% destination states according to wind intensity

% Create waterfall map 
waterMap = zeros(nRow,nCol);
waterMap(:,[3 4]) = 2;
waterMap(:,[5 6]) = 1;
waterRegionIdx = find(waterMap);

for actCounter = 1:length(ActionsList)
    action = ActionsList(actCounter);
    stateTransition = envWaterDeterministic.T(:,:,action2idx(envWaterDeterministic,action)); 
    
    %zoom on affected region by two-subscript indexing
    zoomStateTrans = stateTransition(:,waterRegionIdx);
    %find possible transition in water region
    [currentStateIdx,nextDefaultStateIdx] = find(zoomStateTrans);    
    %apply bias to indexes from zooming
%     currentStateIdx = currentStateIdx + waterRegionIdx(1) - 1; % linear index of current state
    nextDefaultStateIdx = nextDefaultStateIdx + waterRegionIdx(1) - 1; % linear index of next state
    
    %clear default transitions from any start state in water zone
    stateTransition(currentStateIdx,nextDefaultStateIdx) = 0;
    
    %create new transitions for water zone
    [nextDefaultRow,nextDefaultCol] = ind2sub([nRow nCol],nextDefaultStateIdx); 
    nextStateRow = min(nextDefaultRow + waterMap(nextDefaultStateIdx),nRow); 
    nextNewStateIdx = sub2ind([nRow nCol],nextStateRow,nextDefaultCol); 
    newTransitionIdx = sub2ind([stateSize stateSize],currentStateIdx,nextNewStateIdx);  
    stateTransition(newTransitionIdx) = 1;
    
    %apply rule to state transition    
    envWaterDeterministic.T(:,:,action2idx(envWaterDeterministic,action)) = stateTransition;
end

%% Set up obstacles and incorporate to state transitions
% envWaterDeterministic.ObstacleStates = "[3,4]";
updateStateTranstionForObstacles(envWaterDeterministic);

env = rl.env.rlMDPEnv(envWaterDeterministic);


    
    
            

    
