classdef AgentWrapper < matlab.System & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.SampleTime & ...
        matlab.system.mixin.Nondirect
% AGENTWRAPPER: System object used as a wrapper for RL Agents

% Copyright 2017-2018 The MathWorks, Inc.
    properties (Nontunable)
        % root level model
        
        % ModelName
        ModelName = 'model_name'
        
        % the agent will be resolved during setup via the object name
        % (which should reside in the mask workspace of the block that owns
        % this system object)
        
        % Agent Object Name
        AgentObjectName = ''; 
        
        % the block path of the system object needed to resolve the agent
        % object
        
        % Wrapper Block Path
        BlockPath = '';
        
        % treat the system as direct feedthrough. In most cases this should
        % be true, but in situations with model references it can be
        % difficult to manage direct feedthrough
        
        % Treat as direct feed through
        TreatAsDirectFeedthrough = true
        
        % Is reward signal provided
        IsRewardProvided = true
        
        % Is isdone signal provided
        IsIsDoneProvided = true
        
        % Collect logged signals on step
        CollectLoggedSignals = true
        
        % Collect model states on step
        CollectModelStates = true
    end
    properties (Access = private)
        % Agent object, used for learning
        Agent
        % numer of times the block is executed during a single simulation
        StepCount 
        % action dimensions
        ActionDims = []
        % action data type
        ActionDataType
        % action bus name
        ActionBusName
        % action bus element names (used to order the outgoing bus
        % elements)
        ActionBusElementNames
        % action bus element data types
        ActionBusElementDataTypes
        % struct for output bus actions
        ActionBusStruct
        % observation bus element names
        ObservationBusElementNames
        % cell for observation bus
        ObservationBusCell
    end
    methods 
        function this = AgentWrapper(varargin)
            setProperties(this,nargin,varargin{:});
        end
    end
    methods (Access = protected)
        function varargout = outputImpl(this,varargin)
            % overload to get custom behavior
            
            nextobs = varargin{1};
            reward  = varargin{2};
            isdone  = varargin{3};
            
            % if manually requested to stop training stop the simulation
            if this.Agent.TerminateSimulation
                set_param(this.ModelName,'SimulationCommand','stop');
            end
            
            % set logged sigs to []
            loggedSignals = [];
            
            % check if we need to convert the observations from bus
            % (struct) to cell
            if isstruct(nextobs)
                c = this.ObservationBusCell;
                elnames = this.ObservationBusElementNames;
                for i = 1:numel(elnames)
                    c{i} = nextobs.(elnames{i});
                end
                nextobs = c;
            end
            
            if this.StepCount > 0
                % call step on the agent
                action = step(this.Agent,nextobs,reward,isdone,loggedSignals);
            else
                % call getInitialAction on the agent if this is the first time
                % the block is being executed. This will properly
                % initialize the LastAction and LastObservation
                % properties on the agent
                action = getInitialAction(this.Agent,nextobs);
            end
            
            % check if the output is a bus or not
            if isempty(this.ActionBusStruct)
                % cast to the appropriate data type
                varargout{1} = cast(action,this.ActionDataType);
            else
                if ~iscell(action)
                    action = {action};
                end
                % put the actions into the bus structure
                s = this.ActionBusStruct;
                elnames = this.ActionBusElementNames;
                dts = this.ActionBusElementDataTypes;
                for i = 1:numel(elnames)
                    f = elnames{i};
                    s.(f) = cast(action{i},dts{i});
                end
                varargout{1} = s;
            end
            
            % stop the simulation if isdone is true
            if isdone > 0
                set_param(this.ModelName,'SimulationCommand','stop');
            end
        end
        function updateImpl(this,varargin)
            % iterate the step count
            this.StepCount = this.StepCount + 1;
        end
        function setupImpl(this)
            
            % the rl agent block is only supported in normal and
            % accelerator mode
            mode = get_param(this.ModelName,'SimulationMode');
            if ~ismember(mode,{'normal','accelerator'})
                error(message('rl:block:RLAgentInvalidSimMode'));
            end

            % get the agent via slResolve
            this.Agent = slResolve(this.AgentObjectName,this.BlockPath);
            
            if ~isempty(this.Agent)
                % get the action info from the agent
                ainfo = getActionInfo(this.Agent);
                % check if the action info is sourced from a bus
                [isbus,this.ActionBusName] = isBusSpec(ainfo);
                if isbus
                    this.ActionDims = [1 1];
                    this.ActionBusStruct = Simulink.Bus.createMATLABStruct(this.ActionBusName);
                    this.ActionDataType = this.ActionBusName;
                    this.ActionBusElementNames = cellstr([ainfo.Name]');
                    this.ActionBusElementDataTypes = cellstr([ainfo.DataType]');
                    if numel(fields(this.ActionBusStruct)) ~= numel(ainfo)
                        error(message('rl:block:RLAgentBusActionMismatch'));
                    end
                else
                    if ~isscalar(ainfo)
                        error(message('rl:block:RLAgentActionInfoMustBeScalarForNonBus'));
                    end
                    this.ActionDims = ainfo.Dimension;
                    this.ActionBusStruct = [];
                    this.ActionDataType = ainfo.DataType;
                    this.ActionBusElementNames = {};
                end

                % get the observation info
                oinfo = getObservationInfo(this.Agent);
                % check if the observation info is sourced from a bus
                [isbus,obsBusName] = isBusSpec(oinfo); %#ok<ASGLU>
                               
                if isbus
                    this.ObservationBusElementNames = cellstr([oinfo.Name]');
                    this.ObservationBusCell = cell(1,numel(this.ObservationBusElementNames));
                else
                    if ~isscalar(oinfo)
                        error(message('rl:block:RLAgentObservationInfoMustBeScalarForNonBus'));
                    end
                    this.ObservationBusElementNames = {};
                    this.ObservationBusCell = {};
                end
            end
            
            % initialize any states
            reset(this);
        end
        function resetImpl(this)
            % Initialize / reset discrete-state properties
            this.StepCount = 0;
        end
        function checkForSetup(this)
            % check to see if the block needs to be re-setup
            if isempty(this.Agent) || isempty(this.ActionDims)
                setupImpl(this);
            end
        end
        function releaseImpl(~)
        end
        function num = getNumInputsImpl(~)
            % Define total number of inputs for system with optional inputs
            num = 3;
        end
        function num = getNumOutputsImpl(~)
            % Define total number of outputs for system with optional
            % outputs
            num = 1;
        end
        function loadObjectImpl(this,s,wasLocked)
            % Needed for fast restart
            
            % Set public properties and states
            loadObjectImpl@matlab.System(this,s,wasLocked);
            % call setup after load to make sure the action struct is
            % properly initialized
            checkForSetup(this);
        end

        function s = saveObjectImpl(this)
            % Needed for fast restart
            
            % Set public properties and states
            s = saveObjectImpl@matlab.System(this);
        end

        function varargout = getInputNamesImpl(~)
            % Return input port names for System block
            varargout{1} = 'observation';
            varargout{2} = 'reward';
            varargout{3} = 'isdone';
        end

        function varargout = getOutputNamesImpl(~)
            % Return output port names for System block
            varargout{1} = 'action';
        end
        function varargout = isInputDirectFeedthroughImpl(this,varargin)
            % when TreatAsDirectFeedthrough is false, implicit unit delays
            % are added to the input signals. In most cases
            % TreatAsDirectFeedthrough should be true
            
            nu = nargin(this);
            varargout = cell(1,nu);
            for i = 1:nu
                varargout{i} = this.TreatAsDirectFeedthrough;
            end
        end
        function sts = getSampleTimeImpl(this)
            % get the block sample time
            checkForSetup(this);
            ts = getSampleTime(this.Agent);
            if ~isnumeric(ts) || ~isscalar(ts) || ~isfinite(ts) || ~isreal(ts) || ts <= 0
                if isequal(ts,-1)
                    error(message('rl:block:RLAgentSampleTimeCannotBeInherited'));
                else
                    error(message('rl:block:RLAgentInvalidSampleTime'));
                end
            end
            sts = createSampleTime(this,'Type','Discrete','SampleTime',ts);
        end
        function varargout = getOutputSizeImpl(this)
            % get the size and sample time from the agent
            checkForSetup(this);
            varargout{1} = this.ActionDims;
        end
        function varargout = isOutputComplexImpl(this)
            ny = nargout(this);
            varargout = cell(1,ny);
            for i = 1:ny
                varargout{i} = false;
            end
        end
        function varargout = getOutputDataTypeImpl(this)
            checkForSetup(this);
            varargout{1} = this.ActionDataType;
        end
        function varargout = isOutputFixedSizeImpl(this)
            ny = nargout(this);
            varargout = cell(1,ny);
            for i = 1:ny
                varargout{i} = true;
            end
        end
    end
    methods(Access = protected,Static)
        function simMode = getSimulateUsingImpl
            % interpreted execution is only allowed
            simMode = "Interpreted execution";
        end
    end
end