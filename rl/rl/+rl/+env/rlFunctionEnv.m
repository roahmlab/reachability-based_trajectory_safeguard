classdef rlFunctionEnv < rl.env.MATLABEnvironment
% rlFunctionEnv: Create a MATLAB RL environment by supplying step and reset
% functions
%
%   ENV = rlFunctionEnv(OINFO,AINFO,STEPFUN,RESETFUN) 
%   creates an RL environment with specified observation specification
%   OINFO, action specification AINFO, custom step and reset functions 
%   (STEPFUN and RESETFUN). STEPFUN and RESETFUN can be function names, 
%   function handles, or anonymous functions.
%
%   STEPFUN specifies how the environment advances to the next state from a
%   given action. STEPFUN requires the following signature:
%   [Observation, Reward, IsDone, LoggedSignals] = myStepFunction(Action, LoggedSignals)
%
%   RESETFUN sets the default state of the environment. RESETFUN 
%   requires the following signature:
%   [InitialObservation, LoggedSignals] = myResetFunction()
%
%   To pass information from one step to the next, such as the environment 
%   state, use LoggedSignals.
%
%   STEPFUN and RESETFUN must be on the MATLAB path, in the current working
%   folder, or local functions.
%
%   Examples
%   Assume that the custom step and reset functions, 'myStepFunction' and 
%   'myResetFunction', are already defined.
%       % Specify observation and action information
%       ObservationInfo = rlNumericSpec([4 1]);
%       ActionInfo = rlFiniteSetSpec([1 2]);
%
%       % Create an environment using a function handle for the reset 
%       % function and an anonymous function for the step function.
%       StepFun = @(action,loggedSignals) myStepFunction(action,loggedSignals);
%       ResetFun = @myResetFunction;
%       env = rlFunctionEnv(ObservationInfo,ActionInfo,StepFun,ResetFun);
%
%       % Create an environemnt by specifying the function names.
%       env = rlFunctionEnv(ObservationInfo,ActionInfo,'myStepFunction','myResetFunction');
%
%   See also: rlCreateEnvTemplate

    
    % Copyright 2018-2019 The MathWorks Inc.
    
    %% Properties
    properties
        % Specify and initialize environment's necessary properties
        StepFcn
        ResetFcn
    end
    
    properties (SetAccess = private, SetObservable = true)
        % Store user-defined signals (e.g. state)
        LoggedSignals
    end
    
    properties (Access = private)
        % Flag indicates constructor is called
        IsConstructed = false
    end
    
    %% Necessary Methods
    methods
        % Contructor
        function this = rlFunctionEnv(ObservationInfo,ActionInfo,StepFcn,ResetFcn)
            narginchk(4,4)
            this = this@rl.env.MATLABEnvironment(ObservationInfo,ActionInfo);
            % register and validate custom functions
            this.ResetFcn = ResetFcn;
            this.StepFcn = StepFcn;
            this.IsConstructed = true;
            % validate env
            validateEnvironment(this);
        end
        
        function [Observation, Reward, IsDone, LoggedSignals] = step(this, Action)
            % Execute custom step function to advance reinforcement 
            % learning environment based on given action
            %   [Observation, Reward, IsDone, LoggedSignals] = step(Action, LoggedSignals) 
            %   advances the environnment to the next state based on 
            %   action Action. To pass information from one step to the 
            %   next, such as the environment state, use LoggedSignals.

            [Observation, Reward, IsDone, LoggedSignals] = feval(this.StepFcn,Action,this.LoggedSignals);
            this.LoggedSignals = LoggedSignals;
        end
        
        
        function InitialObservation = reset(this)
            % Execute custom reset function to set reinforcement 
            % learning environment to the default state
            %   [InitialObservation, LoggedSignals] = reset()
            %   sets the environnment to the default state. To pass 
            %   information from one step to the next, such as the 
            %   environment state, use LoggedSignals.
            
            InitialObservation = feval(this.ResetFcn);
            %this.LoggedSignals = LoggedSignals;
        end
    end
    
    %% Get/Set Methods
    methods
        % StepFcn
        function set.StepFcn(this,Value)
            % Register custom step function to the environment
            
            % validate StepFcn attribute
            validateattributes(Value,{'char','string','function_handle'},{},'','StepFcn');
            % validate function name/handle, cannot validate anonymous function
            try
                StepNargin = nargin(Value);
            catch ex
                me = MException(message('rl:env:errFuncEnvStepCannotSet'));
                me = addCause(me,ex);
                throw(me);
            end
            if StepNargin ~= 2
                error(message('rl:env:customStepInputArg'));
            end
            this.StepFcn = Value;
            % validate ObservationInfo by simulating environment, only
            % after construction
            if this.IsConstructed
                validateEnvironment(this);
            end
        end
        % ResetFcn
        function set.ResetFcn(this,Value)
            % Register custom reset function to the environment
            
            % validate ResetFcn attribute
            validateattributes(Value,{'char','string','function_handle'},{},'','ResetFcn');
            % validate function name/handle, cannot validate anonymous function
            try
                ResetNargin = logical(nargin(Value));
            catch ex
                me = MException(message('rl:env:errFuncEnvResetCannotSet'));
                me = addCause(me,ex);
                throw(me);
            end
            if ResetNargin
                error(message('rl:env:customResetInputArg'));
            end
            this.ResetFcn = Value;
            % validate ObservationInfo by simulating environment, only
            % after construction
            if this.IsConstructed
                validateEnvironment(this);
            end
        end
    end
end