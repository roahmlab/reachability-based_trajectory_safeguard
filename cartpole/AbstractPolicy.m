classdef AbstractPolicy < matlab.mixin.Copyable
% ABSTRACTPOLICY

% Copyright 2018 The MathWorks, Inc.

    properties (Hidden,SetAccess = protected,SetObservable)
        % flag used to determine if simulation/training should end early
        TerminateSimulation (1,1) logical = false
    end
    properties (Hidden)
        MaxSteps (1,1) {mustBePositive} = Inf
    end
    properties (Access = protected)
        % unit delays to create experiences
        LastObservation = []
        LastAction      = []
        % episodic info
        EpisodeInfo (1,1) rl.util.EpisodeInfo = rl.util.EpisodeInfo();
    end
    properties (Access = private)
        % logger for gathering experiences
        Logger = []
    end
    methods (Hidden)
        function [lastObservation,lastAction] = getLastData(this)
            lastObservation = this.LastObservation;
            lastAction = this.LastAction;
        end
        function terminateSimulation(this,val)
            % request to stop the simulation
            this.TerminateSimulation = val;
            % flush the event queue
            drawnow();
        end
        function info = getEpisodeInfo(this)
            % get episode info from agent
            info = this.EpisodeInfo;
        end
        function exp = checkForEarlyTermination(this,exp)
            % check if the experience was due to the max number of steps
            % being reached
            if ~isinf(this.MaxSteps) && (exp{5} == 0) && (this.EpisodeInfo.StepsTaken == this.MaxSteps)
                exp{5} = 2;
            end
        end
        %% logger functions
        function attachLogger(this,maxSteps)
            % attach a logger 
            this.Logger = rl.util.ExperienceLogger(maxSteps);
        end
        function val = isLoggerAttached(this)
            val = ~isempty(this.Logger) && isvalid(this.Logger);
        end
        function detatchLogger(this)
            % detatch the logger, but don't delete
            this.Logger = [];
        end
        function resetLogger(this)
            if isLoggerAttached(this)
                reset(this.Logger);
            end
        end
        function experiences = getExperiences(this)
            % get the logged experiences
            if isLoggerAttached(this)
                experiences = getExperiences(this.Logger);
            else
                experiences = [];
            end
        end
        %% pre/post sim functions
        function preSimFcn(this,simCount) %#ok<INUSD>
            % overload to do work before simulation starts
            
            % make sure the logger is reset at the top of the simulation
            resetLogger(this);
            % initialize episode info
            reset(this.EpisodeInfo);
        end
        function postSimFcn(this,simCount,simInfo) %#ok<INUSD>
            % overload to do work after simulation ends
        end
    end
    methods
        function action = step(this,observation,reward,isdone,loggedSignals,varargin)
            % ACTION = step(AGENT,OBSERVATION,REWARD,ISDONE) advances the
            % AGENT and returns ACTION from OBSERVATION, REWARD and ISDONE 
            % terminal flag. During training, the AGENT update its 
            % parameters. During simulation, the AGENT returns optimal 
            % ACTION from its policy.
            try
                if isempty(this.LastAction)
                    % set the last action from getActionImpl
                    this.LastAction = getActionImpl(this,localCellify(observation));
                    this.LastAction = iUnwrapDlarray(this.LastAction);
                end
                if isempty(this.LastObservation)
                    % set the last action from the observation
                    this.LastObservation = observation;
                end
                
                % update the episode "states"
                update(this.EpisodeInfo,reward);
                
                % step the policy and create the experience
                action = this.LastAction;
                
                exp = localCellifyExperience({this.LastObservation,action,reward,observation,isdone});
                exp = checkForEarlyTermination(this,exp);

                this.LastObservation = observation;
                % call the step implementation with the generated experience
                action = stepImpl(this,exp);
                action = iUnwrapDlarray(action);
            catch ex  
                me = MException(message('rl:agent:errInvalidExperience'));
                me = addCause(me,ex);
                throw(me);
            end
            this.LastAction = action;
            
            % attach the experience to the logger if available
            if isLoggerAttached(this)
                addExperience2Buffer(this.Logger,{exp});
            end
        end
        function action = getAction(this,observation)
            % ACTION = GETAACTION(THIS,OBSERVATION)
            try
                action = getActionImpl(this,localCellify(observation));
                action = iUnwrapDlarray(action);
            catch ex
                me = MException(message('rl:agent:errInvalidObservation'));
                me = addCause(me,ex);
                throw(me);
            end
        end
        function action = getInitialAction(this,observation)
            % ACTION = GETINITIALACTION(THIS,OBSERVATION)
            % getInitialAction should be called at the top of simulation to
            % initialize the LastObservation and LastAction states. 
            try
                this.LastObservation = observation;
                action = getInitialActionImpl(this,localCellify(observation));
                action = iUnwrapDlarray(action);
                this.LastAction = action;
            catch ex
                me = MException(message('rl:agent:errInvalidObservation'));
                me = addCause(me,ex);
                throw(me);
            end
        end
        function reset(this)
            % RESET(AGENT) resets the AGENT before training.
            this.LastObservation = [];
            this.LastAction      = [];
            % call the reset implementation
            resetImpl(this);
        end
    end
    methods (Abstract)
        ts = getSampleTime(this)
        info = getActionInfo(this)
        info = getObservationInfo(this)
    end
    methods (Abstract,Access = protected)
        % Methods that must be defined by concrete subclasses
        action = stepImpl(this,experience)
        action = getActionImpl(this,observation)
        action = getInitialActionImpl(this,observation)
        resetImpl(this)
    end
    methods (Access = protected)
        function that = copyElement(this)
            that = copyElement@matlab.mixin.Copyable(this);
            that.EpisodeInfo = copy(this.EpisodeInfo);
            if ~isempty(this.Logger)
                that.Logger = copy(this.Logger);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Local Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function exp = localCellifyExperience(exp)
    % observations and actions will be stored as a cell array
    exp{1} = localCellify(exp{1});
    exp{2} = localCellify(exp{2});
    exp{4} = localCellify(exp{4});
    % make sure isdone is single
    exp{5} = single(exp{5});
end
function val = localCellify(val)
    if iscell(val)
        % force to be a row vector cell
        val = val(:)';
    else
        val = {val};
    end
end
function action = iUnwrapDlarray(action)
if iscell(action)
    DlarrayCellIdx = cellfun(@(x) isa(x,'dlarray'), action);
    action(DlarrayCellIdx) = cellfun(@(x) extractdata(x),action(DlarrayCellIdx),'UniformOutput',false);
else
    if isa(action,'dlarray')
        action = extractdata(action);
    end
end
end