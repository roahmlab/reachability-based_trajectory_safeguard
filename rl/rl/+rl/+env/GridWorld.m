classdef GridWorld < rl.env.AbstractMDP
    %GridWorld Creates a grid world
    %  env = rl.env.GridWorld(m,n)
    %
    
    %   Copyright 2018-2019 The MathWorks, Inc.
    
    properties (SetAccess = immutable)
        
        % GridSize
        % [m,n] vector where m and n are the number of rows and columns of
        % the grid world, respectively.
        GridSize
        
    end
    
    %% Dependent Properties
    properties (Dependent)
        
        % CurrentState
        % String with name of current state
        CurrentState string
        
        % States
        % String vector of the state names
        States string
        
        % Actions
        % String vector with the action names
        Actions string
        
        % T
        % State Transisition Matrix T(s,s',a) = P(s'|s,a)
        T
        
        % R
        % Reward Transisition Matrix  r = R(s,s',a)
        R
        
        % ObstacleStates
        % String vector with the state names which can not be reached in
        % the grid world.
        ObstacleStates string
        
        % TerminalStates
        % String vector with the state names which are terminal in the grid
        % world.
        TerminalStates string
     
    end
    
    %% Protected Properties
    properties (Access = protected)
        
        % Obstacle states of the MDP
        ObstacleStates_ = string.empty(0,1)
        
    end
    
    %% Hidden Properties
    properties (Hidden=true)
        
        % Border
        % True if edges moves are not allowed
        Border = true;
        
        % Handle to the GridWorldVisualizer
        Viewer
               
    end
    
    
    %% Public Methods
    methods
        
        function obj = GridWorld(m,n,Moves)
            % GridWorld Construct an m by n grid world with m rows and n
            % columns
            narginchk(2,3)
            
            if ~isnumeric(m) || ~isnumeric(n) || (rem(m,1) ~= 0) || (rem(n,1) ~= 0)
                % REVISIT
                error(message('rl:env:GWValidGridDims'))
            end
            
            if (nargin) == 2
                % Define action space
                Actions = ["N";"S";"E";"W"];
            else
                if strcmpi(Moves,'Standard')
                    Actions = ["N";"S";"E";"W"];
                elseif strcmpi(Moves,'Kings')
                    Actions = ["N";"S";"E";"W";"NE";"NW";"SE";"SW"];
                else
                    error(message('rl:env:GWValidMovesInput'))
                end
            end
            % Define state names
            StateNames = rl.env.GridWorld.createStateNames(m,n);
            % Call super-constructor
            obj = obj@rl.env.AbstractMDP(StateNames,Actions);
            % Cache grid size
            obj.GridSize = [m,n];
            % Initialize state transition matrix
            initializeStateTransitions(obj)
        end
        
        function T = getStateTransition(obj,Actions)
            % Get state transisition matrix
            % T = getStateTransition(obj)
            % T = getStateTransition(obj,Action)
            if nargin == 1
                T = getStateTransition_(obj);
            else
                T = getStateTransition_(obj,Actions);
            end
        end
              
        function RT = getRewardTransition(obj,Action)
            if nargin == 1
                RT = getRewardTransition_(obj);
            else
                RT = getRewardTransition_(obj,Action);
            end
        end
        
        function [S,R,isTerm,loggedsig] = move(obj,Action)
            S0 = obj.CurrentState;
            if isValidAction(obj,S0,Action)
                [S,R,isTerm] = move_(obj,Action);
            else
                error(message('rl:env:MDPInvalidAction'));
            end
            loggedsig = [];
        end
        
        function showViewer(obj)
            % showViewer(env) opens the GridWorld Viewer
            obj.Viewer = rl.env.GridWorldVisualizer(obj.GridSize,[]);
            setTerminalStates(obj.Viewer,obj.TerminalStates)
            setObstacles(obj.Viewer,obj.ObstacleStates)
            updateViewer(obj)
            addlistener(obj,'CurrentStateChanged',@(es,ed) updateViewer(obj));
        end
        
        function updateViewer(obj)
            CS = obj.CurrentState_;
            CS = regexprep(regexprep(CS,"("',"["),")","]");
            obj.Viewer.setCurrentState(eval(CS));
        end
        
        function updateStateTranstionForObstacles(obj)
            ST = obj.StateTransitions_;
            m = obj.GridSize(1);
            n = obj.GridSize(2);
            
            for ct = 1:numel(obj.ObstacleStates)
                % Obstacle
                idxO = state2idx(obj,obj.ObstacleStates(ct));
                
                % N Action
                idxAct = 1;
                idxS1 = find(ST(:,idxO,idxAct));
                for ct1 = 1:numel(idxS1)
                    if idxO+1 < m*n
                        ST(idxS1(ct1),idxO+1,idxAct) = ...
                            ST(idxS1(ct1),idxO+1,idxAct) + ST(idxS1(ct1),idxO,idxAct);
                        ST(idxS1(ct1),idxO,idxAct) = 0;
                    end
                end
                
                % S Action
                idxAct = 2;
                idxS1 = find(ST(:,idxO,idxAct));
                for ct1 = 1:numel(idxS1)
                    if idxO-1 > 0
                        ST(idxS1(ct1),idxO-1,idxAct) = ...
                            ST(idxS1(ct1),idxO-1,idxAct) + ST(idxS1(ct1),idxO,idxAct);
                        ST(idxS1(ct1),idxO,idxAct) = 0;
                    end
                end
                
                
                % E Action
                idxAct = 3;
                idxS1 = find(ST(:,idxO,idxAct));
                for ct1 = 1:numel(idxS1)
                    if idxO-m > 0
                        ST(idxS1(ct1),idxO-m,idxAct) = ...
                            ST(idxS1(ct1),idxO-m,idxAct) + ST(idxS1(ct1),idxO,idxAct);
                        ST(idxS1(ct1),idxO,idxAct) = 0;
                    end
                end
                
                % W Action
                idxAct = 4;
                idxS1 = find(ST(:,idxO,idxAct));
                for ct1 = 1:numel(idxS1)
                    if idxO+m < m*n
                        ST(idxS1(ct1),idxO+m,idxAct) = ...
                            ST(idxS1(ct1),idxO+m,idxAct) + ST(idxS1(ct1),idxO,idxAct);
                        ST(idxS1(ct1),idxO,idxAct) = 0;
                    end
                end
            end
            
            obj.StateTransitions_=ST;
        end
        

        function State = reset(obj,State)
            if nargin == 1
                State = getRandomInitialState(obj);
            else
                PotentialStates = setdiff(obj.States,obj.ObstacleStates);
                PotentialStates = setdiff(PotentialStates,obj.TerminalStates);
                if ~ismember(State, PotentialStates)
                    % REVISIT do we want to change wording of error message
                    error(message('rl:env:MDPInvalidState'))
                end
            end
            obj.CurrentState = State;
        end
    end
    
    
    %% Set/Get Methods
    methods
        % CurrentState
        function set.CurrentState(obj,State)
            % Set Current State
            
            % Validate that its a valid state and not an obstacle state
            obj.CurrentState_ = State;
        end
        
        function State = get.CurrentState(obj)
            % Get Current State
            State = obj.CurrentState_;
        end
        
        % States
        function set.States(obj,States)
            % Set States
            % REVISIT
            error(message('rl:env:GWCannotSetStates'))
        end
        
        function States = get.States(obj)
            % Get States
            States = obj.States_;
        end
        
        % Actions
        function set.Actions(obj,Actions)
            % Set Actions
            % REVISIT
            error(message('rl:env:GWCannotSetActions'))
        end
        
        function Actions = get.Actions(obj)
            % Get Actions
            Actions = obj.Actions_;
        end
        
        % T
        function T = get.T(obj)
            % Get T
            T = getStateTransition(obj);
        end
        
        function set.T(obj,T)
            % Set T
            validateT(obj,T)
            setStateTransition_(obj,T);
        end
        
        %R
        function R = get.R(obj)
            % Get R
            R = getRewardTransition(obj);
        end
        
        function set.R(obj,R)
            % Set R
            validateR(obj,R)
            setRewardTransition_(obj,R);
        end
       
        % ObstacleStates
        function set.ObstacleStates(obj,StateNames)
            for ct = 1:numel(StateNames)
                if isTerminalState(obj,StateNames(ct)) || ~isState(obj,StateNames(ct))
                    % REVISIT
                    error(message('rl:env:GWInvalidObstacleState'))
                end
            end
            obj.ObstacleStates_ = unique(StateNames);
        end
        
        function ObstacleStates = get.ObstacleStates(obj)
            ObstacleStates = obj.ObstacleStates_;
        end
        
        % TerminalStates
        function set.TerminalStates(obj,StateNames)
            for ct = 1:numel(StateNames)
                if isObstacleState(obj,StateNames(ct)) || ~isState(obj,StateNames(ct))
                    % REVISIT
                    error(message('rl:env:GWInvalidTerminalState'))
                end
            end
            obj.TerminalStates_ = unique(StateNames);
        end
        
        function TeriminalState = get.TerminalStates(obj)
            TeriminalState = obj.TerminalStates_;
        end
        
    end
    
    
    %% Protected Methods
    methods (Access = protected)
        
        function initializeStateTransitions(obj)
            for ct = 1:numel(obj.Actions_)
                ActionName = obj.Actions_(ct);
                ST = createStateTransition(obj,ActionName);
                setStateTransition_(obj,ST,ActionName)
            end
        end
        
        function ST = createStateTransition(obj,Action)
            m = obj.GridSize(1);
            n = obj.GridSize(2);
            ST = zeros(m*n,m*n);
            switch Action
                case "N"
                    % Move North
                    for c = 1:n
                        r = 1;
                        if obj.Border
                            ST(r+(c-1)*m,r+(c-1)*m) = 1;
                        else
                            ST(r+(c-1)*m,c*m) = 1;
                        end
                        for r = 2:m
                            ST(r+(c-1)*m,(r-1)+(c-1)*m) = 1;
                        end
                    end
                case "S"
                    % Move South
                    for c = 1:n
                        for r = 1:m-1
                            ST(r+(c-1)*m,(r+1)+(c-1)*m) = 1;
                        end
                        r = m;
                        if obj.Border
                            ST(r+(c-1)*m,r+(c-1)*m) = 1;
                        else
                            ST(c*m,1+(c-1)*m) = 1;
                        end
                    end
                case "E"
                    % Move East
                    for r = 1:m
                        for c = 1:n-1
                            ST(r+(c-1)*m,r+(c+0)*m) = 1;
                        end
                        c = n;
                        if obj.Border
                            ST(r+(c-1)*m,r+(c-1)*m) = 1;
                        else
                            ST(r+(c-1)*m,r) = 1;
                        end
                    end
                case "W"
                    % Move West
                    for r = 1:m
                        c = 1;
                        if obj.Border
                            ST(r,r) = 1;
                        else
                            ST(r,r+(n-1)*m) = 1;
                        end
                        for c = 2:n
                            ST(r+(c-1)*m,r+(c-2)*m) = 1;
                        end
                    end
                case "NW"
                    % Move North West
                    for c = 2:n
                        for r = 2:m
                            ST(r+(c-1)*m,(r-1)+(c-2)*m) = 1;
                        end
                    end
                    if obj.Border
                        for c = 1:n %N
                            r = 1;
                            ST(r+(c-1)*m,r+(c-1)*m) = 1; 
                        end
                        for r = 1:m %W
                            ST(r,r) = 1;                 
                        end
                    end
                case "SW"
                    % Move South West
                    for c = 2:n
                        for r = 1:m-1
                            ST(r+(c-1)*m,(r+1)+(c-2)*m) = 1;
                        end
                    end
                    if obj.Border
                        for c = 1:n %S
                            r = m;
                            ST(r+(c-1)*m,r+(c-1)*m) = 1;
                        end
                        for r = 1:m %W
                            ST(r,r) = 1;
                        end
                    end
                case "NE"
                    % Move North East
                    for c = 1:n-1
                        for r = 2:m
                            ST(r+(c-1)*m,(r-1)+(c+0)*m) = 1;
                        end
                    end
                    if obj.Border
                        for c = 1:n %N
                            r = 1;
                            ST(r+(c-1)*m,r+(c-1)*m) = 1; 
                        end
                        for r = 1:m %E
                            c = n;
                            ST(r+(c-1)*m,r+(c-1)*m) = 1;
                        end
                    end
                case "SE"
                    % Move South East
                    for c = 1:n-1
                        for r = 1:m-1
                            ST(r+(c-1)*m,(r+1)+(c+0)*m) = 1;
                        end
                    end
                    if obj.Border
                        for c = 1:n %S
                            r = m;
                            ST(r+(c-1)*m,r+(c-1)*m) = 1;
                        end
                        for r = 1:m %E
                            c = n;
                            ST(r+(c-1)*m,r+(c-1)*m) = 1;
                        end
                    end
            end
        end
        
        function b = isObstacleState(obj,State)
            % Check if State is an ObstacleState
            b = any(strcmp(State, obj.ObstacleStates));
        end
        
        function State = getRandomInitialState(obj)
            % get random initial state that is not an obstacle or terminal
            % state
            PotentialStates = setdiff(setdiff(obj.States,obj.ObstacleStates),obj.TerminalStates);
            % Randomly select one of the potential states
            State = PotentialStates(randi(numel(PotentialStates)));
            
        end
        
    end
    
    
    %% Static Methods
    methods (Static = true)
        
        function Names = createStateNames(m,n)
            % Create array of state name "[m,n]"   
            
            Names = "";
            ct = 1;
            for ct1 = 1:n
                for ct2 = 1:m
                    Names(ct,1) = sprintf("[%d,%d]",ct2,ct1);
                    ct = ct+1;
                end
            end
        end
        
    end
end

