function env = WaterFallStochasticGridWorld 
% WaterFallStochasticGridWorld: Creates Stochastic Water Fall RL environment.

% Copyright 2018 The MathWorks, Inc.

%% Set up basic GridWorld environment
nRow = 8;    % number of rows
nCol = 7;   % number of columns
envWaterStochastic = createGridWorld(nRow,nCol);%,'Kings'); % REVISIT
startState = "[5,1]";
positiveTermStates = "[4,5]";
negativeTermStates = ["[8,1]";"[8,2]";"[8,3]";"[8,4]";"[8,5]";"[8,6]";"[8,7]"];
terminalStates = [positiveTermStates; negativeTermStates];
envWaterStochastic.CurrentState = startState;
envWaterStochastic.TerminalStates = terminalStates;
ActionsList = envWaterStochastic.Actions;        % get list of possible action
StatesList = envWaterStochastic.States;          % get list of states
nS = numel(StatesList);    % number of states
nA = numel(ActionsList);   % number of possible actions 

%% Incorporate notion of border and obstacles to state transitions
updateStateTranstionForObstacles(envWaterStochastic); %REVISIT: is this needed as it is called later

%% Set Rewards Transitions
envWaterStochastic.R = -1*ones(nS,nS,nA); 
envWaterStochastic.R(:,state2idx(envWaterStochastic,positiveTermStates),:) = 10;
envWaterStochastic.R(:,state2idx(envWaterStochastic,negativeTermStates),:) = -10;

%% Set up waterfall rule
% Create waterfall map with wind intesity as specified in example 6.5
% From default gridworld, find transitions in waterfall zone and change
% destination states according to wind intensity

% Create waterfall map 
waterMap = zeros(nRow,nCol);
waterMap(:,[3 4]) = 2;
waterMap(:,[5 6]) = 1;
waterRegionIdx = find(waterMap);

for i = 1:length(ActionsList)
    action = ActionsList(i);
    stateTransition = envWaterStochastic.T(:,:,action2idx(envWaterStochastic,action));
        
    %zoom on affected region by two-subscript indexing
    zoomStateTrans = stateTransition(:,waterRegionIdx);
    %find possible transition in water region
    [currentStateIdx,nextDefaultStateIdx] = find(zoomStateTrans);   
    %apply bias to indexes from zooming
    nextDefaultStateIdx = nextDefaultStateIdx + waterRegionIdx(1) - 1; % linear index of next state
    % If in waterfall zone, new next state shift up 1 row
    %clear default transitions from any start state in water zone
    stateTransition(currentStateIdx,nextDefaultStateIdx) = 0;
    [nextInitRow,nextInitCol] = ind2sub([nRow nCol],nextDefaultStateIdx); % Find subscripts in gridworld 7x10
    
    % apply stochastic to ST
    nextStateBaseRow = min(nextInitRow + waterMap(nextDefaultStateIdx),nRow); % row if move accordingly to deterministic water value
    nextStateAboveRow = max(nextStateBaseRow-1,1); % result in one row above
    nextStateBelowRow = min(nextStateBaseRow+1,nRow); % result in one row below
    nextStateCol = nextInitCol; 
    % Set 1/3 chance to land on expected state from wind influence
    nextBaseState = sub2ind([nRow nCol],nextStateBaseRow,nextStateCol);
    nextBaseStateIdx = sub2ind([nS nS],currentStateIdx,nextBaseState);
    stateTransition(nextBaseStateIdx) = stateTransition(nextBaseStateIdx) + 1/3;
    % Set 1/3 chance to land 1 row up from expected state from wind influence
    nextAboveState = sub2ind([nRow nCol],nextStateAboveRow,nextStateCol);
    nextAboveStateIdx = sub2ind([nS nS],currentStateIdx,nextAboveState);
    stateTransition(nextAboveStateIdx) = stateTransition(nextAboveStateIdx) + 1/3;
    % Set 1/3 chance to land 1 row down from expected state from wind influence
    nextBelowState = sub2ind([nRow nCol],nextStateBelowRow,nextStateCol);
    nextBelowStateIdx = sub2ind([nS nS],currentStateIdx,nextBelowState);
    stateTransition(nextBelowStateIdx) = stateTransition(nextBelowStateIdx) + 1/3;
    
    envWaterStochastic.T(:,:,action2idx(envWaterStochastic,action)) = stateTransition;
end
%% Set up obstacles and incorporate to state transitions
updateStateTranstionForObstacles(envWaterStochastic);

env = rl.env.rlMDPEnv(envWaterStochastic);


    
    
            

    
