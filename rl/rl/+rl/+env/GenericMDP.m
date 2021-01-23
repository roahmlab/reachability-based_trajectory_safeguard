classdef GenericMDP < rl.env.AbstractMDP
    %GenericMDP Creates a MDP model
    %  env = rl.env.GenericMDP(NumStates,NumActions)
    %  env = rl.env.GenericMDP(StateNames,ActionName)
    %
    %
    
    %   Copyright 2018-2019 The MathWorks, Inc.
       
    %% Dependent Properties
    % REVISIT should all of these be dependent properties?
    properties (Dependent)
        
        % CurrentState
        % String with name of current state
        CurrentState 
        
        % States
        % String vector of the state names
        States 
        
        % Actions
        % String vector with the action names
        Actions 
        
        % T
        % State Transisition Matrix T(s,s',a) = P(s'|s,a)
        T
        
        % R
        % Reward Transisition Matrix  r = R(s,s',a)
        R
        
        % TerminalStates
        % String vector with the state names which are terminal in the grid
        % world.
        TerminalStates 
     
    end
    
    
    %% Hidden Properties
    properties (Hidden=true)
              
        % Handle to the MDP Visualizer
        Viewer
               
    end
    
    
    %% Public Methods
    methods
        
        function obj = GenericMDP(States,Actions)
            % GenericMDP Construct an generic MDP
            %
            % MDP = rl.env.GenericMDP(3,2) 
            %     creates an MDP with 3 states and 2 actions with default
            %     state and action names
            %
            % MDP = rl.env.GenericMDP(["state1";"state2";"state3"],["action1","action2"]) 
            %     creates an MDP with 3 states and 2 actions with the
            %     specified state and action names
            narginchk(2,2)
            
            % States
            if isnumeric(States) 
                if (rem(States,1) == 0)
                    StateNames = string(rl.env.AbstractMDP.strseq('s',1:States));
                else
                    error(message('rl:env:MDPEnvInvalidInput'))
                end
            elseif isstring(States) && isvector(States) && isequal(numel(unique(States)),numel(States))
                StateNames = States(:);
            else
                error(message('rl:env:MDPEnvInvalidInput'))
            end
            
            % Actions
            if isnumeric(Actions) 
                if (rem(Actions,1) == 0)
                    ActionNames = string(rl.env.AbstractMDP.strseq('a',1:Actions));
                else
                    error(message('rl:env:MDPEnvInvalidInput'))
                end
            elseif isstring(Actions) && isvector(Actions) && isequal(numel(unique(Actions)),numel(Actions))
                ActionNames = Actions(:);
            else
                error(message('rl:env:MDPEnvInvalidInput'))
            end

            
            % Call super-constructor
            obj = obj@rl.env.AbstractMDP(StateNames,ActionNames);
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
        
        function State = reset(obj,State)
            if nargin == 1
                State = getRandomInitialState(obj);
            else
                PotentialStates = setdiff(obj.States,obj.TerminalStates);
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
            State = localConvertChar2Str(State);
            if isState(obj,State)
                % Validate state
                obj.CurrentState_ = State;
            else
                error(message('rl:env:MDPInvalidState'))
            end
        end
        
        function State = get.CurrentState(obj)
            % Get Current State
            State = obj.CurrentState_;
        end
        
        % States
        function set.States(obj,States)
            % Set States
            States = localConvertChar2Str(States);
            if isstring(States) && isvector(States) && isequal(numel(unique(States)),numel(obj.States_))
                idx = strcmp(obj.CurrentState_,obj.States_);
                obj.States_ = States(:);
                % Update current state based on new state names
                obj.CurrentState_ = obj.States_(idx);
            else
                error(message('rl:env:MDPEnvInvalidStateActionInput',numel(obj.States_)))
            end
        end
        
        function States = get.States(obj)
            % Get States
            States = obj.States_;
        end
        
        % Actions
        function set.Actions(obj,Actions)
            % Set Actions
            Actions = localConvertChar2Str(Actions);
            if isstring(Actions) && isvector(Actions) && isequal(numel(unique(Actions)),numel(obj.Actions_))
                obj.Actions_ = Actions(:);
            else
                error(message('rl:env:MDPEnvInvalidStateActionInput',numel(obj.Actions_)))
            end
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
        
        % R
        function R = get.R(obj)
            % Get R
            R = getRewardTransition(obj);
        end
        
        function set.R(obj,R)
            % Set R
            validateR(obj,R)
            setRewardTransition_(obj,R);
        end
        
        % TerminalStates
        function set.TerminalStates(obj,StateNames)
            StateNames = localConvertChar2Str(StateNames);
            for ct = 1:numel(StateNames)
                if ~isState(obj,StateNames(ct))
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
    
        function State = getRandomInitialState(obj)
            % get random initial state that is not an obstacle or terminal
            % state
            PotentialStates = setdiff(obj.States,obj.TerminalStates);
            % Randomly select one of the potential states
            State = PotentialStates(randi(numel(PotentialStates)));
            
        end
        
    end
    
    
    %% Static Methods
    methods (Static = true)
       
    end
end

function value = localConvertChar2Str(value)
    if ischar(value) || iscellstr(value) %#ok<ISCLSTR>
        value = string(value);
    end
end

