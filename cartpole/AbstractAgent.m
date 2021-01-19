classdef (Abstract) AbstractAgent < rl.agent.AbstractPolicy
    %RLAGENTABSTRACT  Create abstract class of RL agent.    
    
    % Copyright 2018 The MathWorks, Inc.
    
    %% PROPERTIES
    properties (Access = protected)
        % Observation information data for Agent
        ObservationInfo
        
        % Action information data for Agent
        ActionInfo
        
        % Local copy of training options
        TrainingOptions (1,1) rl.option.rlTrainingOptions
        
        % Does agent have a critic
        HasCritic = false;
        
        % Sample time of the agent
        SampleTime = -1;
    end
    properties (Access = private)
        % stepping mode
        StepMode string = "sim"
        
        % step function handle
        StepFcn function_handle
        
        % intial action function handle
        InitialActionFcn function_handle
        
        % pre sim function handle
        PreSimFcn 
        
        % post sim function handle
        PostSimFcn 
    end
    properties (Hidden)
        % Custom Reward function function handle that returns reward
        % r = CustomRewardFcn(exp)
        CustomRewardFcn
    end
    
    %% ABSTRACT METHODS  
    methods (Abstract, Access = public)
        % compute action with exploration
        action = getActionWithExploration(this,observation);
    end
    methods (Abstract, Access = protected)
        % learn from current experiences, return action with exploration
        % exp = {state,action,reward,nextstate,isdone}
        action = learn(this,exp);
        % validate training options and agent compatibility
        trainOptions = validateAgentTrainingCompatibilityImpl(this,trainOptions)
        % set/get tunable parameters
        p = getLearnableParametersImpl(this)
        setLearnableParametersImpl(this,p)
        % indicate whether the agent uses RNN
        HasState = hasStateImpl(this)
    end
    %% CONCRETE METHODS
    methods
        % constructor
        function this = AbstractAgent()
            % initilize the training options property with the default so
            % we can call sim from Simulink
            this.TrainingOptions = rlTrainingOptions();
            % initialize the step and init action function handles
            setStepMode(this,"sim");
        end          
        %% Get & set methods
        function Info = getObservationInfo(this)
            Info = this.ObservationInfo;
        end
        function Info = getActionInfo(this)
            Info = this.ActionInfo;
        end
        function Ts = getSampleTime(this)
            Ts = this.SampleTime;
        end
        %% codegen
        function varargout = generatePolicyFunction(this,varargin)
            % GENERATEPOLICYFUNCTION(NAME1,VALUE1,...) Generate a
            % standalone m-function that evaluates the learned policy of
            % the agent. The generated function can be used to generate
            % code with MATLAB Coder or GPU Coder to deploy the learned
            % policy. The following name-value pairs are supported:
            %
            %   FunctionName    Name of the generated function. 
            %                   (default = 'evaluatePolicy')
            %   PolicyName      Name of the policy variable in the generated function.
            %                   (default = 'policy')
            %   MATFileName     Name of the generated MAT-file that the generated function will load from.
            %                   (default = 'agentData.mat')
            
            % codegen for RNN agent is not supported
            if hasState(this)
                error(message('rl:agent:errPolicyCodegenStatefulNotSupport'))
            end
            
            p = inputParser;
            addParameter(p,'FunctionName','evaluatePolicy',...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','FunctionName'));
            addParameter(p,'PolicyName','policy',...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','PolicyName'));
            addParameter(p,'MATFileName','agentData.mat',...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','MATFileName'));

            parse(p,varargin{:});
            argStruct = p.Results;
            
            % use Name propery of action/observation info to determine names
            na = numel(this.ActionInfo);
            no = numel(this.ObservationInfo);
            
            ustr = sprintf('observation%u,',1:no);
            ystr = sprintf('action%u,'     ,1:na);
            ustr(end) = '';
            ystr(end) = '';
            
            argStruct.InputString  = ustr;
            argStruct.OutputString = ystr;
            argStruct.EvaluateInputString  = ustr;
            argStruct.EvaluateOutputString = ystr;
            
            % by default, call the evaluate function
            argStruct.EvaluateFunctionName = 'localEvaluate';
            argStruct.LocalFunctionString = '';
            argStruct.BodyString = sprintf('%s = localEvaluate(%s);',...
                rl.codegen.handleMultipleOutputStrings(ystr),ustr);
            
            % get the representation from the concrete agent. Also let the
            % agent modify argStruct to modify pre/post process
            [rep,argStruct] = generateProcessFunctions_(this,argStruct);
            
            % create the args
            fs = fields(argStruct);
            N = 2*numel(fs);
            args = cell(1,N);
            vals = struct2cell(argStruct);
            
            args(1:2:(N-1)) = fs  ;
            args(2:2: N   ) = vals;
            
            % generate the code from the representation
            try
                [varargout{1:nargout}] = generatePolicyFunction(rep,args{:});
            catch ex
                throwAsCaller(ex);
            end
        end
    end
    methods (Access = protected)
        % policy methods
        action = stepImpl(this,experience)
        action = getActionImpl(this,observation)
        
        % helpers
        function TargetRepresentation = updateTarget(this,Representation,TargetRepresentation,SmoothingFactor,TargetUpdateFrequency)
            % update target representation with current representation
            
            % targetUpdateFrequency > 1 yields "periodic" updates
            % smoothingFactor < 1 yields "soft" updates
            % targetUpdateFrequency > 1 & smoothingFactor < 1 yields "periodic-soft" updates
            if mod(this.EpisodeInfo.StepsTaken,TargetUpdateFrequency) == 0
                TargetRepresentation = syncParameters(TargetRepresentation,Representation,SmoothingFactor);
            end
        end
        
        function q0 = evaluateQ0(this,exp) %#ok<INUSD>
            % overload for agents that implement a critic
            q0 = 0;
        end
        function action = getInitialActionImpl(this,observation)
            % get the action at the top of simulation
            
            % If stepFcn || InitialActionFcn are empty, set to sim mode
            % before any simulations starts. This is only necessary to
            % support simulink simulations via the play button after agent
            % load.
            if isempty(this.InitialActionFcn) || isempty(this.StepFcn)
                setStepMode(this,"sim");
            end
            
            action = this.InitialActionFcn(this,observation);
        end
        function [rep,argStruct] = generateProcessFunctions_(this,argStruct) %#ok<INUSD,STOUT>
            % this method must be overloaded for code generation support
            error(message('rl:general:AgentCodeGenNotSupported',class(this)));
        end
        function TrainOptions = validateAgentTrainingCompatibility(this,TrainOptions)
            % agent-training compatibility check
            % not support parallel training for agent that has state
            if ~strcmpi(TrainOptions.Parallelization,'none') && hasState(this)
                error(message('rl:general:errParallelTrainStatefulNotSupport'))
            end
            TrainOptions = validateAgentTrainingCompatibilityImpl(this,TrainOptions);
        end
        function HasState = hasState(this)
            % indicate whether the agent has state
            HasState = hasStateImpl(this);
        end
        %% step and initial action methods
        function action = simStepFcn(this,exp)
            action = getAction(this,exp{4});
        end
        function action = simWithExplorationStepFcn(this,exp)
            action = getActionWithExploration(this,exp{4});
        end
        function action = simWithExplorationAndAppendStepFcn(this,exp)
            action = getActionWithExploration(this,exp{4});
            appendExperience(this,exp);
        end
        function action = learnStepFcn(this,exp)
            action = learn(this,exp);
        end
        function action = simInitActionFcn(this,observation)
            action = getAction(this,observation);
        end
        function action = simWithExplorationInitActionFcn(this,observation)
            action = getActionWithExploration(this,observation);
        end
        function that = copyElement(this)
            that = copyElement@rl.agent.AbstractPolicy(this);
        end
    end
    methods (Access = protected, Static)
        function validateActorCriticInfo(Rep1, Rep2)
            % error check for data spec consistency between 2
            % representations and whether they are both have state or
            % stateless
            % REVISIT: should this be a method on representation?
            
            if ~isempty(Rep1) && ~isempty(Rep2)
                % check for data spec consistency 
                if ~isCompatible(Rep1.ObservationInfo, Rep2.ObservationInfo)
                    error(message('rl:agent:errObservationInfoAC'))
                end
                if (isprop(Rep1,'ActionInfo') && ~isempty(Rep1.ActionInfo)) && ...
                        (isprop(Rep2,'ActionInfo') && ~isempty(Rep2.ActionInfo))
                    if ~isCompatible(Rep1.ActionInfo, Rep2.ActionInfo)
                        error(message('rl:agent:errActionInfoAC'))
                    end
                end
                
                % if stateful training, both representation must be stateful
                if (hasState(Rep1)  && ~hasState(Rep2)) || ...
                        (~hasState(Rep1) &&  hasState(Rep2))
                    error(message('rl:agent:errRepresentationNotAllStateful'))
                end
            end
        end
    end
    methods (Hidden)
        %% PRE/POST TRAIN FOR TRAININGMANAGER
        function preSettings = preTrain(this)
            % overload for custom behavior before training starts
            preSettings = struct();
            % reset the agent
            reset(this);
        end
        function postTrain(this,preSettings) %#ok<INUSD>
            % overload for custom behavior after training ends
        end
        %% pre/post sim functions
        function preSimFcn(this,simCount)
            preSimFcn@rl.agent.AbstractPolicy(this,simCount);
            if ~isempty(this.PreSimFcn)
                this.PreSimFcn(this,simCount);
            end
        end
        function postSimFcn(this,simCount,simInfo) 
            postSimFcn@rl.agent.AbstractPolicy(this,simCount,simInfo);
            if ~isempty(this.PostSimFcn)
                this.PostSimFcn(this,simCount,simInfo);
            end
        end
        %% HELPERS
        function setStepFcn(this,fcn)
            validateattributes(fcn,"function_handle","scalar","setStepFcn","fcn");
            if nargin(fcn) ~= 2
                error(message('rl:agent:errInvalidStepFcn'));
            end
            this.StepMode = "custom";
            this.StepFcn = fcn;
        end
        function setInitialActionFcn(this,fcn)
            validateattributes(fcn,"function_handle","scalar","setInitialActionFcn","fcn");
            if nargin(fcn) ~= 2
                error(message('rl:agent:errInvalidInitialActionFcn'));
            end
            this.StepMode = "custom";
            this.InitialActionFcn = fcn;
        end
        function setPreSimFcn(this,fcn)
            if ~isempty(fcn)
                validateattributes(fcn,"function_handle","scalar","setPreSimFcn","fcn");
                if nargin(fcn) ~= 2
                    error(message('rl:agent:errInvalidPreSimFcn'));
                end
            end
            this.PreSimFcn = fcn;
        end
        function setPostSimFcn(this,fcn)
            if ~isempty(fcn)
                validateattributes(fcn,"function_handle","scalar","setPostSimFcn","fcn");
                if nargin(fcn) ~= 3
                    error(message('rl:agent:errInvalidPostSimFcn'));
                end
            end
            this.PostSimFcn = fcn;
        end
        function setStepMode(this,mode)
            % put the agent in the specified mode
            validatestring(mode,["sim","sim-with-exploration","learn","sim-with-exploration-and-append"],"setStepMode","StepMode");
            switch mode
                case "sim"
                    % compute action with no noise
                    setStepFcn(this,@(agent,exp) simStepFcn(agent,exp));
                    setInitialActionFcn(this,@(agent,obs) simInitActionFcn(agent,obs));
                case "sim-with-exploration"
                    % compute action with noise which is useful when running on a
                    % worker
                    setStepFcn(this,@(agent,exp) simWithExplorationStepFcn(agent,exp));
                    setInitialActionFcn(this,@(agent,obs) simWithExplorationInitActionFcn(agent,obs));
                case "sim-with-exploration-and-append"
                    % compute action with noise and append to experience buffer
                    % useful for parallel worker
                    setStepFcn(this,@(agent,exp) simWithExplorationAndAppendStepFcn(agent,exp));
                    setInitialActionFcn(this,@(agent,obs) simWithExplorationInitActionFcn(agent,obs));
                case "learn"
                    % compute action and learn from experiences
                    setStepFcn(this,@(agent,exp) learnStepFcn(agent,exp));
                    setInitialActionFcn(this,@(agent,obs) simWithExplorationInitActionFcn(agent,obs));
                otherwise
                    % should never be hit
                    error(message('rl:agent:errInvalidStepMode',mode));
            end
            this.StepMode = mode;
        end
        function mode = getStepMode(this)
            mode = this.StepMode;
        end
        function val = hasCritic(this)
            val = this.HasCritic;
        end
        function learnFromExperiences(this,experiences)
            % can be called outside of agent
            
            for i = 1:size(experiences,1)
                exp = experiences{i};
                % learn
                learn(this,exp);
            end
        end
        function experienceBuffer = preprocessExperience(this,experienceBuffer) %#ok<INUSL>
            % Compute representation loss variables (e.g. advantage, critic
            % targets, etc.), and batch inputs from experience buffer
        end
        function gradientBuffer = accumulateGradient(this,bufferTime) %#ok<INUSD>
            % accumulate gradient from experience buffer
            gradientBuffer = [];
        end
        function appendExperience(this,experience) %#ok<INUSD>
            % append experience to experience buffer (if any)
        end
        function applyGradient(this,gradientBuffer) %#ok<INUSD>
            % update representation from gradients
        end
        function p = getLearnableParameters(this)
            % get tunable parameters
            p = getLearnableParametersImpl(this);
        end
        function setLearnableParameters(this,p)
            % set tunable parameters
            setLearnableParametersImpl(this,p)
        end
        function trainMgr = buildTrainingManager(this,env,trainingOptions)
            % attach the training options to the agent and construct a
            % training manager.
            this.TrainingOptions = trainingOptions;
            trainMgr = rl.train.TrainingManager(env,this,trainingOptions);
        end
    end
end