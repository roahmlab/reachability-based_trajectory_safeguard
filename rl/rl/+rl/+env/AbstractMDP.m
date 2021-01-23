classdef AbstractMDP < handle
    %rl.env.AbstractMDP Abstract class for representing an MDP
    %   This class represents the primary features of the MDP
      
    % Copyright 2018-2019 The MathWorks Inc.
        
    %% Public Properties
    properties
        
    end
    
    %% Protected Properties
    properties (Access = protected)
        % First action is empty represents direct transisitions(e.g. pass
        % throughs)
        Actions_ = string.empty(0,1)
        
        % States of the MDP
        States_ = string.empty(0,1)
        
        % Terminal States of the MDP
        TerminalStates_ = string.empty(0,1)
        
        % Current State of the MDP
        CurrentState_ = string.empty(0,1)
        
        % State transition probablity matrix T(s,s',a) = P(s'|s,a)
        % Matrix is size (s,s',a)
        StateTransitions_
        
        % Reward matrix R(s,s',a) = P(r,s'|s,a)?
        % Matrix is size (s,s',a)
        Rewards_
        
        % Recorded History
        RecordHistory_ = false;
        History_  % (s,a,r+1,s+1)
        
        % REVISIT
        RewardFunction_  % REVISIT: IS THIS NEEDED?
    end
    
    %% GET/SET Methods
    methods
        
        function set.CurrentState_(obj,value)
            %Set method for CurrentState_ property
            obj.CurrentState_ = value;
            % Send notification that current state has changed
            notify(obj,'CurrentStateChanged')
        end
        
    end
    
    
    %% Public Methods
    methods (Access = public)
        
        function StateName = idx2state(obj,idx)
            % Convert an index into a state name
            try
                StateName = obj.States_(idx(:));
            catch
                % Invalid idx specified
                error(message('rl:env:MDPInvalidIndex'))
            end
        end
        
        function idx = state2idx(obj,StateName)
            % Convert a state name into an index
            try
                [~,idx] = ismember(StateName,obj.States_);
                idx = idx(:);
                if any(idx==0)
                    error(message('rl:env:MDPInvalidState'))
                end
            catch
                % Invalid StateName specified
                error(message('rl:env:MDPInvalidState'))
            end
        end
        
        function ActionName = idx2action(obj,idx)
            % Convert an index into a action
            try
                ActionName = obj.Actions_(idx(:));
            catch
                % Invalid idx specified
                error(message('rl:env:MDPInvalidIndex'))
            end
        end
        
        function idx = action2idx(obj,ActionName)
            % Convert a action name into an index
            try
                [~,idx] = ismember(ActionName,obj.Actions_);
                idx = idx(:);
                if any(idx==0)
                    error(message('rl:env:MDPInvalidAction'))
                end
            catch
                % Invalid StateName specified
                error(message('rl:env:MDPInvalidAction'))
            end
        end

    end
    
    
    %% Protected Methods
    methods (Access = protected)
        function obj = AbstractMDP(StateNames,ActionNames)
            %AbstractMDP Construct an instance of this class
            
            % Revisit when inputs are cellstr
            if nargin == 1
                ActionNames = [];
            end
            obj.States_ = StateNames(:);
            obj.Actions_ = ActionNames(:); 
            nStates = numel(obj.States_);
            obj.StateTransitions_= zeros(nStates,nStates,numel(obj.Actions_));
            obj.Rewards_= zeros(nStates,nStates,numel(obj.Actions_));
            obj.CurrentState_ = StateNames(1);
        end
        
        function addAction_(obj,ActionName, StateTransition)
            %addAction Summary of this method goes here
            
             % REVISIT check ActionName and StateTransisition size
             % REVISIT Validate rows of state transitions and values
             % 0<= v <=1
            obj.Actions_(end+1,1) = ActionName;
            obj.StateTransitions_(:,:,end+1) = StateTransition;
        end
        
        function ST = getStateTransition_(obj, ActionName)
            %getStateTransition_ Summary of this method goes here
            if nargin == 1
                ST = obj.StateTransitions_;
            else
                % REVISIT check ActionName size
                ind = action2idx(obj,ActionName);
                ST = obj.StateTransitions_(:,:,ind);
            end
        end
        
        function setStateTransition_(obj, ST, ActionName)
            %setStateTransition_ Summary of this method goes here
            if nargin == 2
                obj.StateTransitions_=ST;
            else      
                ind = action2idx(obj,ActionName);
                obj.StateTransitions_(:,:,ind)=ST;
            end
        end
        
        function RT = getRewardTransition_(obj,ActionName)
            %getRewardTransition Summary of this method goes here
            if nargin == 1
                RT = obj.Rewards_;
            else
                % REVISIT check ActionName size
                ind = action2idx(obj,ActionName);
                RT = obj.Rewards_(:,:,ind);
            end
        end
        
        function setRewardTransition_(obj, RT, ActionName)
            %setRewardTransition Summary of this method goes here
            if nargin == 2
                obj.Rewards_=RT;
            else
                ind = action2idx(obj,ActionName);
                obj.Rewards_(:,:,ind)=RT;
            end
        end
        
        function R = getReward_(obj,S0,A,S1)
           S0_idx = state2idx(obj,S0); 
           S1_idx = state2idx(obj,S1); 
           RT = getRewardTransition_(obj,A);
           RTall = getRewardTransition_(obj,"");
           
           R = RT(S0_idx,S1_idx)+RTall(S0_idx,S1_idx);
        end
        
        function addStates_(obj, StateNames)
            %addStates_ Summary of this method goes here
            %   Detailed explanation goes here
            n = numel(StateNames);
            for ct = 1:n
                obj.States_(end+1,1) = StateNames(ct);
            end
            z = zeros(n,numel(obj.States_));
            obj.StateTransitions_(:,:,:)
        end
        
        function removeStates_(obj, StateNames)
            %removeStates_ Summary of this method goes here
            %   Detailed explanation goes here
        end
        
        function b = isValidAction(obj, S, ActionName)
            ST = getStateTransition(obj,ActionName);
            b = any(strcmp(obj.States_(any(ST,2)),S));
        end
           
        function [S,R,isTerm] = move_(obj,ActionName)
            %move Summary of this method goes here
            R = 0;
            S0 = obj.CurrentState_;
            isTerm = isTerminalState(obj);
            if isTerm
                R = 0;
                S = S0;
            elseif  true; %isValidAction(obj,S0,ActionName) revisit
                ST = getStateTransition(obj,ActionName);
                S0_idx = state2idx(obj,S0 );
                T = ST(S0_idx,:);
                S1_idx = find(cumsum(T)>=rand,1);
                S = obj.States_(S1_idx);% replace with s1idx
                
                RT = getRewardTransition_(obj,ActionName);
                               
                R = RT(S0_idx,S1_idx);
                
                obj.CurrentState_ = S;
                isTerm = isTerminalState(obj);

            else
                % REVISIT Error
            end
            

        end
        
        function S1 = predict(obj,ActionName,S0)
            % predict Predicts the next state using the transition matrix
            %
            % S1 = predict(obj,ActionName)
            %     Predicts the next state given action using the MDP's
            %     current state
            % S1 = predict(obj,ActionName,S0)
            %     Predicts the next state given action and the specified 
            %     current state

            if nargin == 2
                S0 = obj.CurrentState_;
            end
            ST = getStateTransition(obj,ActionName);
            S0_idx = state2idx(obj,S0 );
            T = ST(S0_idx,:);
            S1_idx = find(cumsum(T)>=rand,1);
            S1 = obj.States_(S1_idx);
        end
        
        function b = isState(obj,State)
            % Check if State is a State of the MDP
            b = any(strcmp(State, obj.States_));
        end
        
        function b = isTerminalState(obj,State)
            % Check if current state is terminal
            if nargin == 1
                State = obj.CurrentState_;
            end
            b = false;
            if ~isempty(obj.TerminalStates_) && any(strcmp(obj.TerminalStates_,State))
                b = true;
            end
        end
        
        function clearHistory(obj)
            % Clears history
            obj.History_ = [];
        end
        
        function record(obj,Flag)
            % Sets the Record History flag
            obj.RecordHistory_ = Flag;
        end
        
        function validateR(obj,R)
            if ~(isnumeric(R) && isequal(size(R),size(obj.Rewards_)))
               error(message('rl:env:MDPInvalidR'))
            end
        end
        
        function validateT(obj,T)
            if ~(isnumeric(T) && isequal(size(T),size(obj.StateTransitions_)) && ...
                    all(T>=0 & T<=1,'all') && ...
                    all(sum(T,2)==1 | sum(T,2)==0,'all'))
                error(message('rl:env:MDPInvalidT'))
            end
        end
    end
    
    methods (Static)
        
        function strvec = strseq(str,idx)
            %STRSEQ  Builds a sequence of indexed strings.
            %
            %   STRVEC = STRSEQ(STR,INDICES) returns the string vector STRVEC obtained
            %   by appending the integer values INDICES to the string STR. For example,
            %    	strseq('e',[1 2 4])
            %   returns
            %   	string({'e1';'e2';'e4'})
            %
            n = numel(idx);
            strvec = cell(n,1); % preallocate
            for i=1:n
                strvec{i,1} = sprintf('%s%d',str,idx(i));
            end
            strvec = string(strvec);
        end
        
    end
    
    
    
    %% Events
    events
        CurrentStateChanged
    end
end

