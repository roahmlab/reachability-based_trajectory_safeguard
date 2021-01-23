classdef rlMDPEnv < rl.env.MATLABEnvironment
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

    % Copyright 2017-2018 The MathWorks Inc.  
    
    properties
        Model rl.env.AbstractMDP
        ResetFcn
    end
    
    
    %% Public Methods
    methods
        function obj = rlMDPEnv(MDP)
            %rlMDPEnv(MDP) Construct an MDP environment for reinforcement learning
            %   MDP should be a rl.env.AbstractMDP
            narginchk(1,1)
            if ~(isa(MDP, 'rl.env.AbstractMDP') && isscalar(MDP))
                error(message('rl:env:errRLMDPEnvInvalidInput'))
            end
            ActionInfo = rlFiniteSetSpec(1:numel(MDP.Actions));
            ActionInfo.Name = 'MDP Actions';
            
            ObservationInfo = rlFiniteSetSpec(1:numel(MDP.States));
            ObservationInfo.Name = 'MDP Observations';
            
            % get Observation and Action information from MDP
            obj = obj@rl.env.MATLABEnvironment(ObservationInfo,ActionInfo); 
           
            obj.Model = MDP;

        end
        
         % set plotting for gridworld
        function plot(obj)
            showViewer(obj.Model);
        end
        
        function state = getCurrentState(obj)
            state = obj.Model.CurrentState;
        end
        
    end
    
    %% Implement Abstact Methods
    methods 
        
        % define step function on gridworld
        function [Observation,Reward,isTerminal,LoggedSignals] = step(this,Action)
            LoggedSignals = [];
            Action = idx2action(this.Model,Action);
            [Observation,Reward,isTerminal] = move(this.Model,Action);
            Observation = state2idx(this.Model,Observation);
        end
        
        % set reset for gridworld
        function State = reset(this,State)
            % REVISIT
            StateSpecified = false;
            if nargin == 1 && ~isempty(this.ResetFcn)
                State = feval(this.ResetFcn);
                StateSpecified = true;
            elseif nargin == 2
                StateSpecified = true;
            end
                
            if StateSpecified
                State = idx2state(this.Model,State);
                State = reset(this.Model,State);
            else
                State = reset(this.Model);
            end
            
            State = state2idx(this.Model,State);

        end

    end
end

