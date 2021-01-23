classdef rlAbstractRepresentation
    % RLREPRESENTATIONBSTRACT: Defines AD interface for RL.
    
    % Copyright 2018 The MathWorks, Inc.
    
    properties
        % Representation Options
        %
        % Specifies options for how the representation is updated, type
        % "help rlRepresentationOptions" for details.
        Options
    end
    
    properties (Access = protected)
        % Description of representation
        rlFunction % REVISIT: is this property needed here? Should this move to classes that need it?
        adModel
        
        % Validate = true if model size and properties been validated
        Validated = false;
        
        % Stores local options for representations
        LocalOptions % REVISIT: is this property needed here? Should this move to classes that need it?
    end

    properties (Hidden) % REVISIT: do these need to be public?
        % Defines the type and size of the observation data spaces 
        %   (array of rl.util.rlDataSpec).
        ObservationInfo %rl.util.RLDataSpec
        
        % Defines the type and size of the action data spaces
        %   (array of rl.util.rlDataSpec).
        ActionInfo %rl.util.RLDataSpec

        % Nameing order for channels of the representation
        ObservationNames
        
        % Nameing order for channels of the representation
        ActionNames     
        
        % Index mapping from ObservationInfo and ObservationNames
        % [Obs1, Obs2, Obs3] is the observatin info. This is the same order
        % as the exp comes in, then if ObservationNames was
        % [Obs2,Obs3,Obs1] then Observation index would be [2,3,1]
        ObservationIndex
        
        % Index mapping from ActionInfo and ActionNames
        ActionIndex
        
        % Maps the combined [obs,act] to the input channels of the
        % representation Q=f(s,a) input [s,a] and inputindex [observationindex, num(obs)+actionindex]
        InputIndex
    end    
    
    %% Constructor
    methods
        % Constructor
        function this = rlAbstractRepresentation(rlFunction,varargin)
            
            this.rlFunction = rlFunction;         
            
            %% set rlRepresentation Options
            ix = find(cellfun(@(x) isa(x,'rl.option.rlRepresentationOptions'),varargin),1);
            if isempty(ix)
                this.Options = rlRepresentationOptions();                
            else
                this.Options = varargin{ix};
            end            
            
            % REVISIT
            if isa(rlFunction,'rlTable')
                return
            end
            
            %% Observation and action name value pairs
            Observation = [];
            ix = find(cellfun(@(x) (ischar(x) || isstring(x)) && strcmpi(x,'Observation'),varargin),1);
            if ~isempty(ix) && ix < length(varargin)
                Observation = varargin{ix+1};                
            end
            if ~isempty(Observation) && (ischar(Observation) || isstring(Observation))
                % convert single char and string to cell array char
                Observation = { char(Observation) };
            end
            if isempty(Observation) || ...
                    ~(iscell(Observation) && all(cellfun(@(x) isstring(x) || ischar(x),Observation)))
                error(message('rl:agent:errObservationNameRequired'));
            end
            this.ObservationNames = reshape(Observation,1,[]);
            
            Action = [];
            ix = find(cellfun(@(x) (ischar(x) || isstring(x)) && strcmpi(x,'Action'),varargin),1);
            if ~isempty(ix) && ix < length(varargin)
                Action = varargin{ix+1};
            end            
            if ~isempty(Action) && (ischar(Action) || isstring(Action))
                % convert single char and string to cell array char
                Action = { char(Action) };
            end            
            if (~isempty(Action) || any(cellfun(@any,(cellfun(@(x) strcmpi(x,'action'),varargin,'UniformOutput',false))))) && ...
                    ~(iscell(Action) && all(cellfun(@(x) isstring(x) || ischar(x),Action)))
                error(message('rl:agent:errActionNameRequired'));
            end
            this.ActionNames = reshape(Action,1,[]);
            
            
            %% Observation and action info
            ix = find(cellfun(@(x) isa(x,'rl.util.RLDataSpec'),varargin),2);
            switch length(ix)
                case 0
                    error(message('rl:agent:errObservationInfoRequired'));
                case 1
                    this.ObservationInfo = reshape(varargin{ix(1)},1,[]);
                    this.ActionInfo = [];
                otherwise
                    this.ObservationInfo = reshape(varargin{ix(1)},1,[]);
                    this.ActionInfo = reshape(varargin{ix(2)},1,[]);
            end
            validateattributes(this.ObservationInfo,{'rl.util.RLDataSpec'},{'vector'},'','ObservationInfo')
            if ~isempty(this.ActionInfo)
                validateattributes(this.ActionInfo,{'rl.util.RLDataSpec'},{'vector'},'','ActionInfo')
            end
            
            %% Action name and info required
            if ~isempty(this.ActionNames) && isempty(this.ActionInfo)
                error(message('rl:agent:errActionInfoRequired'));
            end
            if isempty(this.ActionNames) && ~isempty(this.ActionInfo)
                error(message('rl:agent:errActionNameRequired'));
            end
            
            %% only single channel actions are supported
            if ~this.Options.AllowMultipleActionChannels && ~isempty(this.ActionInfo) && ~isscalar(this.ActionInfo)
                error(message('rl:agent:errMultiActionChannelsNotSupported'));
            end
        end
        function this = set.Options(this,opt)
            if isempty(opt)
                opt = rlRepresentationOptions;
            end

            validateattributes(opt,{'rl.option.rlRepresentationOptions'},{'scalar'},'','Options');
            if ~isequal(this.Options,opt)
                % check GPU available
                if strcmpi(opt.UseDevice,'gpu')
                    try
                        % this function errors if gpu is not available
                        gpuAvailable = nnet.internal.cnn.util.GPUShouldBeUsed(opt.UseDevice); %#ok<NASGU>
                    catch
                        % it's not available, rever to cpu
                        opt.UseDevice = 'cpu';
                        warning(message('rl:agent:errGPUNotAvailable'));
                    end
                end
                
                this.Options = opt;
                % update existing options
                this = updateOptions(this,opt);
            end
        end
    end
    
    %% Abstract Public Methods
    methods (Abstract)
        % Compute gradients of representation output wrt learnable parameters or inputs
        [grad,paramerer] = gradient(this,dy,dx,inputValues,varargin);      
%         function[grad,parameter] = gradient(obj,dy,dx,inputValues,varargin)
%             % This function only supports the following pairs and API
%             % Required: (dOutput/dParameter,dOutput/dAction,dOutput/dObservation,dLoss/dParameter)
%             % dy is 'output' or 'loss'
%             % dx is 'parameter','observation' or 'action'
%             % inputValues is the values passed into the representation
%             % varagin is 
%         end
        
        % Step learnable parameters with gradient information or loss
        this = step(this,gradient,parameter);
        
        % Fit learnable parameters with input and target data
        this = fit(this,Data,Target);        
        
        % Copy representation
        Representation2 = copy(Representation1);                    
   
        % Set the loss function for the network
        this = setLoss(this,AutomaticLoss)

        % MDL = getModel(REP) gets the computation model MDL of function
        % representation REP
        % * Deep Neural Network Representation
        %       MDL = getModel(REP) returns a SeriesNetwork, DAGNetwork or
        %       dlnetwork object.
        % * Table Representation
        %       MDL = getModel(REP) returns an rlTable object.
        % * Linear Basis Representation
        %       MDL = getModel(REP) returns a struct contains the basis
        %       function and the learnable parameters.
        mdl = getModel(this)
    end
    
    %% Abstract Protected Methods
    methods (Access=protected,Abstract)
        Value = evaluate_(this,Data)
        Value = getLearnableParameterValues_(this);
        this = setLearnableParameterValues_(this,Value);
        sz = getSize_(this,identifier);
        validateLearnableParameterSize(this,Value);
    end
    
    %% Public Methods
    methods
        %% Get size of observation, action, output, and parameter
        function sz = getSize(this,identifier)
            if this.Validated
                validstr = validatestring(identifier,{'observation','action','output','input','parameter'});
                sz = getSize_(this,validstr);
            else
                error(message('rl:agent:errLossRequired'));
            end
        end
        
        %% Compute output of representation
        function Value = evaluate(this,Data)
            if this.Validated
                try
                    Value = evaluate_(this,Data);
                catch ME
                    if ~(iscell(Data) && numel(Data)==(length(this.InputIndex)))
                        error(message('rl:agent:errIncompatibleInput1'));
                    end
                    DataSize = cellfun(@(x) size(x),Data,'UniformOutput',false);
                    InputSize = arrayfun(@(x) x.Dimension,this.ObservationInfo,'UniformOutput',false);
                    % add action size only if it appears at input
                    % assumes action does not appear at input and output at the
                    % same time
                    if ~isempty(this.ActionInfo) && strcmp(this.ActionIndex{1},'input')
                        InputSize = [InputSize arrayfun(@(x) x.Dimension,this.ActionInfo,'UniformOutput',false)];
                    end
                    
                    % Check dimensions are consistent
                    if ~all(cellfun(@(x,y) rl.util.isArrayEqual(x,y),InputSize,DataSize))
                        error(message('rl:agent:errIncompatibleInput2'));
                    end
                    
                    % Check batch sizes are consistent
                    if all(cellfun(@(x) length(x)==4,DataSize))
                        if length(unique(cellfun(@(x) x(4),DataSize)))~=1
                            error(message('rl:agent:errIncompatibleInput2'));
                        end
                    end
                    
                    % the error cannot be trapped and rethrown
                    rethrow(ME);                    
                end
            else
                error(message('rl:agent:errLossRequired'));
            end                        
        end        
        
        %% Get/set learnable parameter values
        function Value = getLearnableParameterValues(this)
            % getLearnableParameterValues: RL representation returns the
            % learnable parameter values as a cell array
            %  
            %   val = getLearnableParameterValues(rep) where rep is the
            %   representation.
            
            if this.Validated
                Value = getLearnableParameterValues_(this);
            else
                error(message('rl:agent:errLossRequired'));
            end
        end
        function this = setLearnableParameterValues(this,Value)
            % setLearnableParameterValues: Sets the values of learnable
            % parameters in representation to the values passed as a cell
            % array
            %  
            %   rep = setLearnableParameterValues(rep,val) where rep is the
            %   representation and val is the cell array of values.
            if this.Validated
                validateLearnableParameterSize(this,Value);
                this = setLearnableParameterValues_(this,Value);
            else
                error(message('rl:agent:errLossRequired'));
            end                
        end        
    end
    
    %% Hidden Methods
    methods (Hidden)
        %% Update learnable parameters of one representation with others
        function this = updateParameter(this,Representation,smoothingFactor)
            % validateRepresentationType(Representation) if function is visible
            if nargin<3
                smoothingFactor = 1;
            %else
            %    validateattributes(Value,{'numeric'},{'scalar','real','positive','<=',1},'','smoothingFactor');                
            end            
            ParameterValues = getLearnableParameterValues_(Representation);
            CurrentParameterValues = getLearnableParameterValues_(this);
            for ct=1:numel(CurrentParameterValues)
                CurrentParameterValues{ct} = CurrentParameterValues{ct}*(1-smoothingFactor) + ParameterValues{ct}*smoothingFactor;
            end
            this = setLearnableParameterValues_(this,CurrentParameterValues);
        end        
        function Data = mapInputValues(this,Data)
            [~,ix]=sort(this.InputIndex);
            Data = Data(ix);
        end
        function hasActionInfo(this,Description)
            if isempty(this.ActionInfo) || ~isa(this.ActionInfo,'rl.util.RLDataSpec')
                error(message('rl:agent:errNoActionInfo',sprintf('%s',Description)));
            end
        end
        function validateActionType(this,Description,Type)
            if strcmpi(Type,'discrete')
                errorFlag = ~isa(this.ActionInfo,'rl.util.rlFiniteSetSpec');
            elseif strcmpi(Type,'continuous')
                errorFlag = ~isa(this.ActionInfo,'rl.util.rlNumericSpec');
            end
            if errorFlag
                error(message('rl:agent:errMismatchedActionType',sprintf('%s',Description),sprintf('%s',Type)));
            end
        end
        function validateScalarOutput(this,Description)
            repOutputSize = getSize(this,'output');
            if ~isequal(repOutputSize{1}, [1 1])
                error(message('rl:agent:errRepresentationOutSizeNotScalar',sprintf('%s',Description)))
            end
        end
        function validateRepresentationIOWithAction(this,Description,IO,MatchFlag)
            if MatchFlag
                % positive: representation I/O must be actions
                if ~isempty(this.ActionIndex) && ~strcmpi(this.ActionIndex{1},IO)
                    error(message('rl:agent:errRepresentationActionIOPos',sprintf('%s',IO),...
                        sprintf('%s',Description)));
                end
            else
                % negative: representation I/O cannot be actions
                if ~isempty(this.ActionIndex) && strcmpi(this.ActionIndex{1},IO)
                    error(message('rl:agent:errRepresentationActionIONeg',sprintf('%s',IO),...
                        sprintf('%s',Description)));
                end
            end
        end
        function vals = usampleInputs(this)
            % keep the random seed to restore after usample for
            % reproducability
            s = rng;
            vals = reshape(usample(this.ObservationInfo),1,[]);
            if ~isempty(this.ActionInfo) && strcmpi(this.ActionIndex{1},'input')
                vals = [vals reshape(usample(this.ActionInfo),1,[])];
            end
            rng(s);
        end
        function flag = isActor(this)
            flag = ~isempty(this.ActionInfo) && strcmpi(this.ActionIndex{1},'output');
		end
        function varargout = generatePolicyFunction(this,varargin)
            % generate a stand-alone static policy function. Useful for
            % code generation
            p = inputParser;
            addParameter(p,'PolicyName','policy',...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','PolicyName'));
            addParameter(p,'MATFileName','cgactor.mat',...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','MATFileName'));
            addParameter(p,'FunctionName','cgevaluate',...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','FunctionName'));
            addParameter(p,'InputString','observation',...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','InputString'));
            addParameter(p,'OutputString','action',...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','OutputString'));
            addParameter(p,'EvaluateInputString','input',...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','EvaluateInputString'));
            addParameter(p,'EvaluateOutputString','action',...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','EvaluateOutputString'));
            addParameter(p,'LocalFunctionString','',...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','LocalFunctionString'));
            addParameter(p,'EvaluateFunctionName','localEvaluate',...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','EvaluateFunctionName'));
            addParameter(p,'BodyString','',...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','BodyString'));
            
            parse(p,varargin{:});
            argStruct = p.Results;
            % generate the evaluate function (representation specific)
            argStruct = generateEvaluateFunction_(this,argStruct);
            % generate the code
            [varargout{1:nargout}] = rl.codegen.generateEvaluateFcn(argStruct);
        end
    end
    
    %% Protected Methods
    methods (Access = protected)
        function this = validateRepresentationSizes(this,ObservationSizes,ActionSizes)
            % check observation sizes
            validateRepresentationObservationSizes(this,ObservationSizes)
            % check action sizes
            validateRepresentationActionSizes(this,ActionSizes)
        end
        function validateRepresentationObservationSizes(this,ObservationSizes)
            % get observation info sizes
            ObservationInfoSizes = arrayfun(@(x) x.Dimension,this.ObservationInfo,'UniformOutput',false);
            % check observation sizes
            if numel(ObservationSizes)~=numel(ObservationInfoSizes)
                error(message('rl:agent:errIncompatibleObservation1'));
            end
            if ~all(cellfun(@(x,y) rl.util.isArrayEqual(x,y),ObservationSizes,ObservationInfoSizes))
                error(message('rl:agent:errIncompatibleObservation2'));
            end
        end
        function validateRepresentationActionSizes(this,ActionSizes)
            if ~isempty(ActionSizes)
                ActionInfoSizes = arrayfun(@(x) x.Dimension,this.ActionInfo,'UniformOutput',false);
                if numel(ActionSizes)~=numel(ActionInfoSizes)
                    error(message('rl:agent:errIncompatibleAction1'));
                end
                % check action sizes and action info sizes
                if strcmpi(this.ActionIndex{1},'input')
                    if ~all(cellfun(@(x,y) rl.util.isArrayEqual(x,y),ActionSizes,ActionInfoSizes))
                        error(message('rl:agent:errIncompatibleAction2'));
                    end                                        
                end
                % check the number of action output and a
                if strcmpi(this.ActionIndex{1},'output')
                    ActionInfoDiscreteSizes = [];
                    if isa(this.ActionInfo,'rl.util.rlFiniteSetSpec')
                        ActionInfoDiscreteSizes = arrayfun(@(x) [length(x.Elements) 1],this.ActionInfo,'UniformOutput',false);
                    end
                    if isempty(ActionInfoDiscreteSizes)
                        % continuous action output case
                        if ~all(cellfun(@(x,y) rl.util.isArrayEqual(x,y),ActionSizes,ActionInfoSizes))
                            error(message('rl:agent:errIncompatibleAction2'));
                        end
                    else
                        % discrete action output case
                        if ~all(cellfun(@(x,y) rl.util.isArrayEqual(x,y),ActionSizes,ActionInfoDiscreteSizes))
                            error(message('rl:agent:errIncompatibleAction2'));
                        end
                    end
                end
            end
        end
        function argStruct = generateEvaluateFunction_(this,argStruct) %#ok<INUSD>
            % representations must overload this method for code gen
            % support
            error(message('rl:general:RepCodeGenNotSupported',class(this)));
        end
        function this = updateOptions(this,opt)            
        end
    end
end
