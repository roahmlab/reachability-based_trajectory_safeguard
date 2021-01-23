%rlAbstractRepresentation Defines the interface of reinforcement learning 
% representation. A representation can be value functions, policies or
% a combination of value function and policy.
%
% rlAbstractRepresentation properties:
%    Options         - Representation Option
%    ObservationInfo - Observation data specification
%    ActionInfo      - Action data specification
%    Model           - Computation model object of class rl.representation.model.rlAbstractModel
%
% rlAbstractRepresentation public methods
%    evaluate - Evaluate the model from given input data.
%    gradient - Compute gradient with various types. E.g. dLoss/dParameter
%    setLoss  - Set the model's loss.
%    optimize - Update the model's parameters from an 
%               optimization routine given the gradient values.
%    syncParameters - Update learnable parameters of one representation
%                     with another.
%    get/setLearnableParameters
%    get/setState
%    resetState
%    getSize - Get input/output or parameters size.
%    generatePolicyFunction

% Copyright 2019 The MathWorks, Inc.

classdef rlAbstractRepresentation < matlab.mixin.Heterogeneous
    
    properties (Dependent, SetAccess = private)
        ObservationInfo
    end
    
    properties (Dependent)
        Options
    end
    
    properties (Access = private)
        % Specify options for how the representation is updated, type
        % "help rlRepresentationOptions" for details.
        Options_
    end
    
    properties (Access = protected)
        % Defines the type and size of the observation data spaces, can be
        % an array of rl.util.RLDataSpec
        ObservationInfo_
        
        % Defines the type and size of the action data spaces, can be
        % an array of rl.util.RLDataSpec, can be empty
        ActionInfo_
        
        % Computation representation of the representation
        Model
        
        % Reference for input dimension on each port
        InputDataDimension
        
        % Indicate whether loss is set or not
        IsLossSet
        
        % Solver performs gradient update
        Solver
    end
    
    %======================================================================
    % Public API
    %======================================================================
    methods
        function this = rlAbstractRepresentation(Model,ObservationInfo,ActionInfo,Options)
            % Constructor
            
            this.Model = Model;
            this.ObservationInfo_ = ObservationInfo(:)';
            this.ActionInfo_ = ActionInfo(:)';
            this.Options = Options;
            
            % REVISIT: only single channel actions are supported
            if ~this.Options.AllowMultipleActionChannels && ~isempty(this.ActionInfo_) && ~isscalar(this.ActionInfo_)
                error(message('rl:agent:errMultiActionChannelsNotSupported'));
            end
            
            % cache data dimension to enforce data convention
            % - observation will always be inputs
            % - action dimension will be captured on the representation 
            %   level (only single output Q and determ actor)
            this.InputDataDimension = {this.ObservationInfo.Dimension};
            
            this.IsLossSet = false;
        end
        
        function [Val, State, BatchSize, SequenceLength] = evaluate(this, Data)
            % EVALUATE evaluate the representation with given input data.
            %   [VALUE, STATE] = EVALUATE(REP, DATA) evaluates the
            %   representation REP with input data DATA. It returns output
            %   of the evaluation VALUE and the state of the representation
            %   STATE. 
            %
            %       DATA is a cell array with as many elements as the number
            %       of representation input channels. Each cell array element 
            %       must follow the dimension convention 
            %       SingleDataDimension-by-BatchSize-by-SequenceLength
            %
            %       VALUE is a cell array with as many elements as the number
            %       of representation output channels.
            
            try
                % infer batch size and sequence length to correct dimension
                [BatchSize,SequenceLength] = inferDataDimension(this.ObservationInfo(1), Data{1});
                % evaluate model
                [Val, State] = evaluate(this.Model, Data, BatchSize, SequenceLength);
                % remove singleton dimension of output
                Val = cellfun(@(x) squeeze(x), Val, 'UniformOutput', false);
            catch ex
                validateInputData(this, Data);
                % if the error cannot be trapped and rethrown
                me = MException(message('rl:agent:errRepEvaluate'));
                throw(addCause(me,ex))
            end
        end
        
        function Grad = gradient(this, GradType, InputData, varargin)
            % GRADIENT compute gradients from the representation and input
            % data.
            %   GRAD = GRADIENT(REP, GRADTYPE, INPUTDATA) computes 
            %   gradients from representation REP and input data INPUTDATA.
            %   The gradient type is determined by GRADTYPE.
            %       GRADTYPE can be 'loss-parameters', 'output-input', 'output-parameters'
            %
            %   GRAD = GRADIENT(REP, 'loss-parameters', INPUTDATA, LOSSVARIABLE)
            %   computes the gradients of the representation's loss with
            %   respect to its learnable parameters. Must set loss with
            %   SETLOSS function to use this gradient type.
            %
            %   GRAD = GRADIENT(REP, 'output-input', INPUTDATA)
            %   computes the gradients of the representation's outputs with
            %   respect to its inputs.
            %
            %   GRAD = GRADIENT(REP, 'output-parameters', INPUTDATA)
            %   computes the gradients of the representation's outputs with
            %   respect to its learnable parameters.
            %
            %       INPUTDATA is a cell array with as many elements as the number
            %       of representation input channels. Each cell array element 
            %       must follow the dimension convention 
            %       SingleDataDimension-by-BatchSize-by-SequenceLength
            %
            %       LOSSVARIABLE contains any additional information
            %       required to compute the loss. It is passed as is to
            %       the loss function.
            
            %   GRAD = GRADIENT(REP, 'output-parameters', INPUTDATA, INITIALGRADIENT)
            %   computes the gradients of the representation's outputs from
            %   the given initial gradients with respect to its inputs.
            %
            %       INITIALGRADIENT must have the same data type, shape and
            %       size with the representation output given INPUTDATA.
            
            % REVISIT: add output-action, output-observation? require
            % action/observation input index
            GradType = validatestring(GradType,["loss-parameters","output-input","output-parameters"], 'gradient', 'GradType', 2);
            
            try
                % cast to correct dimension
                [BatchSize,SequenceLength] = inferDataDimension(this.ObservationInfo(1), InputData{1});
                % compute gradient
                Grad = gradient(this.Model, GradType, InputData, BatchSize, SequenceLength, varargin{:});
            catch ex
                if strcmpi(GradType,"loss-parameters") && ~this.IsLossSet
                    error(message('rl:agent:errLossRequiredForGradient'));
                end
                validateInputData(this, InputData);
                % if the error cannot be trapped and rethrown
                if strcmpi(ex.identifier, 'nnet_cnn:internal:cnn:layer:CustomOutputLayer:ForwardLossErrored')
                    % hide trace from internal FcnLossLayer, suggest check
                    % loss function
                    oldCause = ex.cause;
                    ex = MException('',getString(message('rl:agent:errCannotEvaluateLoss')));
                    ex = addCause(ex,oldCause{1});
                end
                me = MException(message('rl:agent:errRepGradient'));
                throw(addCause(me,ex))
            end
        end

        function this = optimize(this, Grad)
            % OPTIMIZE update the learnable parameters of the
            % representation with given gradients following the
            % optimization strategy specified by the representation options
            % at construction.
            %
            %   REP = OPTIMIZE(REP, GRAD) update the learnable parameters 
            %   of the representation REP with gradients GRAD.
            
            LearnableParameters = getLearnableParameters(this);

            % L2 regularize gradients
            if this.Options.L2RegularizationFactor ~= 0
                Grad = rl.internal.optimizer.regularizeGradient(Grad, ...
                    LearnableParameters, this.Options.L2RegularizationFactor);
            end
            
            % threshold gradients
            if ~isinf(this.Options.GradientThreshold)
                Grad = rl.internal.optimizer.thresholdGradient(Grad, ...
                    this.Options.GradientThreshold, ...
                    this.Options.GradientThresholdMethod);
            end
            
            % parameters update
            [this.Solver, LearnableParameters] = calculateUpdate(this.Solver,LearnableParameters,Grad,this.Options.LearnRate);
            
            % set parameters to representation
            this = setLearnableParameters(this,LearnableParameters);
        end
        
        function this = syncParameters(this, Rep, SmoothFactor)
            % SYNCPARAMETERS update the learnable parameters of one
            % representation with another.
            %
            %   TARGETREP = SYNCPARAMETERS(TARGETREP, SOURCEREP, SMOOTHFACTOR) 
            %   updates the learnable parameters of the TARGETREP 
            %   representation with thos of SOURCEREP.
            %   TARGETPARAMS = SMOOTHFACTOR*SOURCEPARAMS + (1 - SMOOTHFACTOR)*TARGETPARAMS
            
            narginchk(2,3)
            if nargin < 3
                SmoothFactor = 1;
            end
            validateattributes(SmoothFactor, {'numeric'}, {'scalar', 'finite', 'nonnan', '<=', 1}, 'syncParameters', 'SmoothFactor', 3);
            
            if SmoothFactor < 1
                this.Model = syncParameters(this.Model, getRLModel(Rep), SmoothFactor);
            else
                % Copy parameters 1-to-1
                this = setLearnableParameters(this, getLearnableParameters(Rep));
            end
        end
        
        function State = getState(this)
            % GETSTATE return the current state of the representation. 
            %
            %   STATE = GETSTATE(REP) returns the current state STATE of 
            %   the representation REP. E.g. current hidden state of a 
            %   recurrent neural network representation.
            
            State = getState(this.Model);
        end
        
        function this = setState(this, State)
            % SETSTATE set the current state of the representation. 
            %
            %   REP = SETSTATE(REP, STATE) sets the current state of the 
            %   representation REP to STATE.
            
            this.Model = setState(this.Model, State);
        end
        
        function this = resetState(this)
            % RESETSTATE reset the state of the representation to zeros.
            %
            %   REP = RESETSTATE(REP) resets the current state of the 
            %   representation REP to zeros.
            
            this.Model = resetState(this.Model);
        end
        
        function this = setLoss(this, Loss)
            %REP = SETLOSS(REP, LOSS) sets loss function handle LOSS to all
            % channel of the representation REP. LOSS can be a character
            % vector, string scalar or a function handle.
            
            % REVISIT: support multiple losses for 1 output
            % REVISIT: support a losses that takes multiple outputs
            
            validateattributes(Loss,{'function_handle','string','char'},{'nonempty'},'setLoss','Loss',1);
            if isa(Loss,'char')
                Loss = string(Loss);
            end
            if isequal(Loss, "")
                error(message('rl:agent:errLossEmptyString'))
            end
            if ~isscalar(Loss)
                error(message('rl:agent:errMultiLoss'))
            end
            if ~isa(Loss,'function_handle')
                Loss = str2func(Loss);
            end
            try
                nargin(Loss);
            catch CauseException
                BaseException = MException(message('rl:agent:errLossCannotSet'));
                BaseException = addCause(BaseException,CauseException);
                throw(BaseException);
            end
            this.Model = setLoss(this.Model, Loss);
            this.IsLossSet = true;
        end
        
        function LearnableParameters = getLearnableParameters(this)
            % GETLEARNABLEPARAMETERS return the learnable parameters of the
            % representation.
            %
            %   LEARNABLEPARAMETERS = GETLEARNABLEPARAMETERS(REP) returns
            %   the learnable parameters LEARNABLEPARAMETERS of the
            %   representation REP.
            
            LearnableParameters = getLearnableParameters(this.Model);
        end
        
        function this = setLearnableParameters(this, LearnableParameters)
            % SETLEARNABLEPARAMETERS set the learnable parameters of the
            % representation.
            %
            %   REP = SETLEARNABLEPARAMETERS(REP, LEARNABLEPARAMETERS) sets
            %   the learnable parameters of the representation REP to
            %   LEARNABLEPARAMETERS value.
            
            this.Model = setLearnableParameters(this.Model, LearnableParameters);
        end
        
        function varargout = generatePolicyFunction(this,varargin)
            % GENERATEPOLICYFUNCTION generate a stand-alone static policy 
            % function. Useful for code generation
            
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
            argStruct = generateEvaluateFunction(this.Model,argStruct);
            % generate the code
            [varargout{1:nargout}] = rl.codegen.generateEvaluateFcn(argStruct);
        end
        
        function Model = getModel(this)
            % MODEL = GETMODEL(REP) returns the computation model MODEL of
            % the representation REP.
            % * Deep Neural Network Model
            %       MODEL = GETMODEL(REP) returns a SeriesNetwork or 
            %       DAGNetwork object.
            % * Table Model
            %       MODEL = GETMODEL(REP) returns an rlTable object.
            % * Basis Function Model
            %       MODEL = GETMODEL(REP) returns a 1x2 cell contains the
            %       basis function handle and the parameters.
            
            Model = getExternalModel(getRLModel(this));
        end
        %==================================================================
        % Get/set
        %==================================================================
        function Options = get.Options(this)
            Options = this.Options_;
        end
        function this = set.Options(this,Options)
            validateattributes(Options,{'rl.option.rlRepresentationOptions'},{'scalar'},'','Option');
            if ~isequal(this.Options,Options)
                % check GPU available
                if strcmpi(Options.UseDevice,'gpu')
                    try
                        % this function errors if gpu is not available
                        gpuAvailable = nnet.internal.cnn.util.GPUShouldBeUsed(Options.UseDevice); %#ok<NASGU>
                        % solver will automatically updates statistics to
                        % gpuArray when interact with model
                    catch
                        % it's not available, revert to cpu
                        Options.UseDevice = 'cpu';
                        warning(message('rl:agent:errGPUNotAvailable'));
                    end
                end
                
                % update model to use cpu or gpu
                if ~isempty(this.Options) && ~strcmpi(this.Options.UseDevice, Options.UseDevice)
                    this.Model = updateExecutionEnvironment(this.Model, Options.UseDevice);
                end
                
                % create Solver if Solver is empty
                if isempty(this.Solver)
                    this.Solver = rl.util.createSolverFactory(Options, this.Model);
                end
                
                % recreate Solver if any of the following options is modified
                % rationale: solver statistics significantly outdated when
                % these options are change and are no longer valid
                if ~isempty(this.Options) && ...
                        (~isequal(this.Options.GradientThreshold, Options.GradientThreshold) ...
                        || ~isequal(this.Options.L2RegularizationFactor, Options.L2RegularizationFactor) ...
                        || ~isequal(this.Options.Optimizer, Options.Optimizer) ...
                        || ~isequal(this.Options.OptimizerParameters, Options.OptimizerParameters))
                    this.Solver = rl.util.createSolverFactory(Options, this.Model);
                end
                
                this.Options_ = Options;
            end
        end
        function ObservationInfo = get.ObservationInfo(this)
            ObservationInfo = this.ObservationInfo_;
        end
    end
    
    %======================================================================
    % Util methods to be called by agent class but not user-facing
    %======================================================================
    methods (Hidden)
        function Model = getRLModel(rep)
            % get RL model of class
            % rl.representation.rlAbstractRepresentation
            
            Model = rep.Model;
        end
        
        function HasState = hasState(this)
            HasState = hasState(this.Model);
        end
        
        function sz = getSize(this,Identifier)
            % Get the representation I/O or parameter size.
            % Use for validation. Return cell array.
            
            ValidStrings = {'input', 'output', 'parameter'};
            Identifier = validatestring(Identifier,ValidStrings);
            sz = getSize(this.Model,Identifier);
        end
    end
    
    %======================================================================
    % Backward compatibility
    %======================================================================
    methods (Hidden)
        function LearnableParameters = getLearnableParameterValues(this)
            % GETLEARNABLEPARAMETERVALUES return the learnable parameters 
            % of the representation.
            %
            %   LEARNABLEPARAMETERS = GETLEARNABLEPARAMETERVALUES(REP) 
            %   returns the learnable parameters LEARNABLEPARAMETERS of the
            %   representation REP.
            
            LearnableParameters = getLearnableParameters(this.Model);
        end
        
        function this = setLearnableParameterValues(this, LearnableParameters)
            % SETLEARNABLEPARAMETERVALUES set the learnable parameters of 
            % the representation.
            %
            %   REP = SETLEARNABLEPARAMETERVALUES(REP, LEARNABLEPARAMETERS)
            %   sets the learnable parameters of the representation REP to
            %   LEARNABLEPARAMETERS value.
            
            this.Model = setLearnableParameters(this.Model, LearnableParameters);
        end
    end
    
    %======================================================================
    % Validation
    %======================================================================
    methods (Access = private)
        function validateInputData(this, Data)
            % Input data error catching for evaluate() and gradient()
            
            % check num data channels == num rep input channels
            if ~(iscell(Data) && numel(Data)==(numel(this.InputDataDimension)))
                error(message('rl:agent:errRepIncorrectNumInputChannel'));
            end
            
            % data dimension on each channel has to follow [DIM BatchSize SequenceLength]
            DataSize = cellfun(@(x) size(x),Data,'UniformOutput',false);
            DataSize = reshape(DataSize, size(this.InputDataDimension));
            
            % the difference of data number of dimension on all channel has
            % to input data spec
            %   0 (BatchSize, SequenceLength = 1)
            %   1 (BatchSize > 1, SequenceLength = 1)
            %   2 (BatchSize > 1, SequenceLength > 1)
            NDimsDiff = cellfun(@(x,y) numel(x) - numel(y), DataSize, this.InputDataDimension);
            if ~all(ismember(NDimsDiff,[0 1 2]))
                error(message('rl:agent:errRepIncorrectInputDim'));
            end
            % check if DIM on each channel is consistent
            if ~all(cellfun(@(x,y) isequal(x(1:numel(y)),y), DataSize, this.InputDataDimension))
                error(message('rl:agent:errRepIncorrectInputDim'));
            end
            
            % check batch, sequence dimension
            if any(NDimsDiff == 2)
                if ~hasState(this)
                    error(message('rl:agent:errNonStatefulRepSequenceData'));
                end
            end
        end
    end
    
    methods (Access = protected)
        function validateModelInputDimension(this)
            % validate model input dimensions with input data dimension
            % spec. e.g. DNN inputImageInput layer size
            % must be called on concrete class constructor
            
            % check if the model has the same number of input channel
            % specified by data specs
            ModelInputSize = getSize(this.Model,'input');
            if numel(ModelInputSize) ~= numel(this.InputDataDimension)
                error(message('rl:agent:errIncompatibleModelInputNumChannel'));
            end
            
            % check if each output channel of the model has compatible size
            % with data specs dimension
            ModelInputSize = reshape(getSize(this.Model,'input'),size(this.InputDataDimension));
            if ~all(cellfun(@(x,y) rl.util.isArrayEqual(x,y), ModelInputSize, this.InputDataDimension))
                error(message('rl:agent:errIncompatibleModelInputDim'));
            end
        end
    end
end