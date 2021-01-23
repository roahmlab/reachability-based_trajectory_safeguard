%rlFunctionModel Implements rlAbstractModel interface. rlFunctionModel
%utilizes autodiff to define function and parameters as RL model.

%   Copyright 2019 The MathWorks, Inc.

classdef rlBasisFunctionModel < rl.representation.model.rlAbstractModel
    
    properties (Access = private)
        % Function handle that generates features from obs
        BasisFcn_
        
        % Parameters unlabelled dlarray
        Parameters_
        
        % Input sizes
        InputSizes_
        
        % Output sizes
        OutputSizes_
        
        % Feature sizes (output of evaluating BasisFcn_)
        FeatureSizes_
        
        % Bias (all zeros) to model linear basis with fullyconnect
        Bias_
        
        % Loss function handle
        LossHandle_
        
        % Output type determines post processing routine apply to the model
        % output. E.g. softmax, softplus
        OutputType_ = "none"
        
        % Output post process function handle, set by sampling strategy
        PostProcessFcn_
    end
    
    methods
        function this = rlBasisFunctionModel(Model,UseDevice,InputSize,OutputSize)
            % BasisModel = rlBasisFunctionModel(Model,InputSize,OutputSize)
            % Output
            %   - BasisModel: rlBasisFunctionModel object, subclass of rlAbstractModel
            % Inputs
            %   - Model     : 2-element cell array of basis function handle
            %                 and intial parameters
            %   - UseDevice : 'cpu' or 'gpu'
            %   - InputSize, OutputSize: cell array of model inputs and
            %                            outputs size.
            
            % REVISIT: support state
            
            narginchk(4,4)
            
            % model must be a 2-element cell array
            if numel(Model) ~= 2
                error(message('rl:agent:errInvalidInputBasisModel'));
            end
            BasisFcn = Model{1};
            Parameters = Model{2};
            if ~isa(BasisFcn,'function_handle') || ~isnumeric(Parameters)
                error(message('rl:agent:errInvalidInputBasisModel'));
            end
            
            % register model sizes
            this.BasisFcn_ = BasisFcn;
            this.InputSizes_ = InputSize;
            
            % get feature size by doing an feval on BasisFcn
            DummyInput = cell(numel(InputSize),1);
            for ct = 1:numel(InputSize)
                DummyInput{ct} = ones(InputSize{ct});
            end
            % basis function must have 1 column vector output
            NumBasisOutputChannel = nargout(this.BasisFcn_);
            if NumBasisOutputChannel > 1
                error(message('rl:agent:errBasisNumOutputArgGt1'))
            end
            Feature = feval(this.BasisFcn_,DummyInput{:});
            if ~iscolumn(Feature)
                error(message('rl:agent:errBasisOutputNotColumnVec'))
            end
            this.FeatureSizes_ = size(Feature);
            NumBasisFeature = size(Feature,1);
            
            % validate output size is a column vector
            this.OutputSizes_  = OutputSize;
            if numel(OutputSize) > 1
                error(message('rl:agent:errBasisFcnModelNumOutputChannelGt1'))
            end
            if ~iscolumn(zeros(OutputSize{1}))
                error(message('rl:agent:errBasisFcnModelOutputNotColumnVec'))
            end
            
            % verify parameter size
            % Out = W'B => size(W) = [nB nOut]
            TrueParamSize = [NumBasisFeature, OutputSize{1}(1)];
            if ~isequal(size(Parameters),TrueParamSize)
                error(message('rl:agent:errBasisFcnModelInvalidParams'))
            end
            this.Parameters_ = dlarray(Parameters);
            this.Bias_ = dlarray(zeros(TrueParamSize(2),1));
            
            this = updateExecutionEnvironmentImpl(this, UseDevice);
        end
    end
    
    %======================================================================
    % Implementation of abstract method
    %======================================================================
    methods (Access = protected)
        function [Val, State] = evaluateImpl(this, InputData, BatchSize, SequenceLength)
            % EVALUATE evaluate the model with given input data.
            %   [VALUE, STATE] = EVALUATE(MODEL, DATA, BATCHSIZE, SEQUENCELENGTH)
            %   evaluates the model MODEL with input data DATA. It returns output
            %   of the evaluation VALUE and the state of the model
            %   STATE.
            %
            %       DATA is a cell array with as many elements as the number
            %       of model input channels. Each cell array element
            %       must follow the dimension convention
            %       SingleDataDimension-by-BatchSize-by-SequenceLength
            %
            %       VALUE is a cell array with as many elements as the number
            %       of model output channels.
            %
            %       BATCHSIZE and SEQUENCELENGTH are scalar given by the
            %       representation level.
            
            if SequenceLength > 1
                error(message('rl:agent:errBasisSequenceData'))
            end
            
            Feature = generateFeatureFromInput(this, InputData, BatchSize);
            % REVISIT: support multiple outputs from BasisFcn_,
            % requires to support multiple action channels
            
            Val = forward(this,Feature,this.Parameters_);
            % return cell to be consistent with interface
            % TODO: determine whether we want to return dlarray or numeric
%             Val = {extractdata(Val)};
            Val = {Val};
            % REVISIT: not support state
            State = {};
        end
        
        function Grad = gradientImpl(this, GradType, InputData, BatchSize, SequenceLength, LossVariable)
            % GRADIENT compute gradients from the model and input
            % data.
            %   GRAD = GRADIENT(MODEL, GRADTYPE, DATA, BATCHSIZE, SEQUENCELENGTH)
            %   computes gradients from model MODEL and input data DATA. The
            %   gradient type is determined by GradType.
            %       GRADTYPE can be 'loss-parameters', 'output-input' or 'output-parameters'
            %
            %   GRAD = GRADIENT(MODEL, 'loss-parameters', DATA, BATCHSIZE, SEQUENCELENGTH, LOSSVARIABLE)
            %   computes the gradients of the model's loss with
            %   respect to its learnable parameters. Must set loss with
            %   SETLOSS function to use this gradient type.
            %
            %   GRAD = GRADIENT(MODEL, 'input-output', DATA, BATCHSIZE, SEQUENCELENGTH)
            %   computes the gradients of the model's outputs with
            %   respect to its inputs.
            %
            %       DATA is a cell array with as many elements as the number
            %       of model input channels. Each cell array element
            %       must follow the dimension convention
            %       SingleDataDimension-by-BatchSize-by-SequenceLength
            %
            %       LOSSVARIABLE contains any additional information
            %       required to computed the loss, usually a struct.
            %       If LOSSVARIABLE is a struct or single cell, it is
            %       passed as is to the loss function.
            %       If LOSSVARIABLE is a cell array with numel == number of
            %       model's outputs, each cell will be passed to the loss
            %       function of each output respectively.
            %
            %       BATCHSIZE and SEQUENCELENGTH are scalar given by the
            %       representation level.
            %
            %   GRAD = GRADIENT(MODEL, 'output-parameters', DATA, BATCHSIZE, SEQUENCELENGTH)
            %   computes the gradients of the model's output with
            %   respect to its learnable parameters.
            
            if SequenceLength > 1
                error(message('rl:agent:errBasisSequenceData'))
            end
            
            switch GradType
                case 'loss-parameters'
                    Feature = generateFeatureFromInput(this, InputData, BatchSize);
                    Grad = {dlfeval(@dLoss_dParams, this, this.Parameters_, Feature, LossVariable)};
                case 'output-input'
                    % cast input data to dlarray to record gradient
                    InputData = iCastDLArray(InputData);
                    Grad = dlfeval(@dOutput_dInput, this, InputData, BatchSize);
                case 'output-parameters'
                    Feature = generateFeatureFromInput(this, InputData, BatchSize);
                    Grad = {dlfeval(@dOutput_dParams, this, this.Parameters_, Feature)};
                otherwise
                    error(message('rl:agent:errInvalidSyntaxGradient'))
            end
        end
        
        function this = setLossImpl(this, Loss)
            %MODEL = SETLOSS(MODEL, LOSS) sets loss function handle LOSS to all
            % channel of the model MODEL. LOSS is a scalar function
            % handle.
            
            this.LossHandle_ = Loss;
        end
        
        function LearnableParameters = getLearnableParametersImpl(this)
            % GETLEARNABLEPARAMETERS return the learnable parameters of the
            % model.
            %
            %   LEARNABLEPARAMETERS = GETLEARNABLEPARAMETERS(MODEL) returns
            %   the learnable parameters LEARNABLEPARAMETERS of the
            %   model MODEL.
            %       LEARNABLEPARAMETERS must satisfy SETLEARNABLEPARAMETERS
            
            LearnableParameters = {this.Parameters_};
        end
        
        function this = setLearnableParametersImpl(this, LearnableParameters)
            % SETLEARNABLEPARAMETERS sets the learnable parameters of the
            % model.
            %
            %   MODEL = SETLEARNABLEPARAMETERS(MODEL, LEARNABLEPARAMETERS) sets
            %   the learnable parameters of the model MODEL to
            %   LEARNABLEPARAMETERS value.
            
            this.Parameters_ = LearnableParameters{1};
        end
        
        function State = getStateImpl(~)
            State = {};
        end
        
        function this = setStateImpl(this, ~)
            % no-op
        end
        
        function this = resetStateImpl(this)
            % no-op
        end
        
        function argStruct = generateEvaluateFunctionImpl(this,argStruct)
            % Generate the evaluate function for the table model
            outputstr = argStruct.EvaluateOutputString;
            inputstr  = argStruct.EvaluateInputString;
            localfcnstr = argStruct.LocalFunctionString;
            evalfcnname = argStruct.EvaluateFunctionName;
            outputstr_y = rl.codegen.handleMultipleOutputStrings(outputstr);
            
            % extract the function
            fcn = this.BasisFcn_;
            wcell = getLearnableParameters(this);
            w = wcell{1};
            fcnstr = func2str(fcn);
            wstr = mat2str(w);
            str1 = sprintf([...
                'fcn = %s;\n',...
                'w = %s;\n'],...
                fcnstr,...
                wstr);
            
            % TODO: continuous stochastic
            switch this.OutputType_
                case "none"
                    str2 = sprintf('%s = w''*fcn(%s);\n',...
                        outputstr_y,inputstr);
                case "probability"
                    % add softmax to the generated code
                    str2 = sprintf('%s = rl.codegen.softmax(w''*fcn(%s));\n',...
                        outputstr_y,inputstr);
                otherwise
                    % this should not happen
                    error(message('rl:agent:errBasisFcnModelOutputTypeNotSupport'))
            end
            
            localevalfcnstr = sprintf([...
                'function %s = %s(%s)\n',...
                '%s%s',...
                'end'],...
                outputstr_y,evalfcnname,inputstr,...
                str1,str2);
            
            % attach local eval fcn string to the arg struct
            if isempty(localfcnstr)
                localfcnstr = localevalfcnstr;
            else
                localfcnstr = sprintf('%s\n%s',localfcnstr,localevalfcnstr);
            end
            argStruct.LocalFunctionString = localfcnstr;
        end
        
        function this = updateExecutionEnvironmentImpl(this,UseDevice)
            % UPDATEEXCUTIONENVIRONMENT updates model to use CPU or GPU
            % based on UseDevice
            
            if strcmpi(UseDevice, 'gpu')
                this.Parameters_ = dlupdate(@gpuArray,this.Parameters_);
                this.Bias_       = dlupdate(@gpuArray,this.Bias_);
            else
                if isa(this.Parameters_,'gpuArray')
                    this.Parameters_ = dlupdate(@gather,this.Parameters_);
                    this.Bias_       = dlupdate(@gather,this.Bias_);
                end
            end
        end
        
        function InternalModel = getInternalModelImpl(this)
            % GETINTERNALMODEL returns the internal model
            
            InternalModel = {this.BasisFcn_ this.Parameters_};
        end
        
        function sz = getSizeImpl(this,Identifier)
            % GETSIZE get the model I/O or parameter size for validation.
            % Identifier can be 'input', 'output', 'parameter'
            
            switch Identifier
                case 'output'
                    sz = this.OutputSizes_;
                case 'input'
                    sz = this.InputSizes_;
                case 'parameter'
                    sz = cellfun(@(x) size(x),getLearnableParameters(this),'UniformOutput',false);
            end
        end
        
        function HasState = hasStateImpl(~)
            % HASSTATE returns true if the model has state (e.g. RNN model)
            
            HasState = false;
        end
        
        function this = setOutputTypeImpl(this,OutputType,varargin)
            % MODEL = SETOUTPUTTYPEIMPL(MODEL, OUTPUTTYPE, POSTPROCESSFCN) 
            % set model output to be compatible with required output type
            %   OUTPUTTYPE: not use for function model
            % e.g. discrete sampling strategy assigns softmax 
            % POSTPROCESSFCN to ensure it outputs probability
            
            this.OutputType_ = OutputType;
            this.PostProcessFcn_ = varargin{1};
        end
        
        function ExternalModel = getExternalModelImpl(this)
            % GETEXTERNALMODEL returns a cell contains the basis function 
            % and the learnable parameters.
            
            ExternalModel = getInternalModel(this);
        end
    end
    
    %======================================================================
    % Gradient methods
    %======================================================================
    methods (Access = private)
        function GradVal = dLoss_dParams(this, Parameters, Feature, LossVariable)
            % dLoss/dParams: gradient of loss wrt parameters
            
            FunctionOutput = forward(this,Feature,Parameters);
            Loss = feval(this.LossHandle_, FunctionOutput, LossVariable);
            GradVal = dlgradient(Loss, Parameters);
        end
        
        function GradVal = dOutput_dInput(this, InputData, BatchSize)
            % dOut/dIn: gradient of output wrt to input
            % REVISIT: add index to choose which input to return
            
            % NOTE: gradient captures feature generation
            Feature = generateFeatureFromInput(this, InputData, BatchSize);
            FunctionOutput = forward(this,Feature,this.Parameters_);
            GradVal = dlgradient(sum(FunctionOutput,'all'), InputData);
        end
        
        function GradVal = dOutput_dParams(this, Parameters, Feature)
            % dOut/dParams: gradient of output wrt to parameters
            
            FunctionOutput = forward(this,Feature,Parameters);
            GradVal = dlgradient(sum(FunctionOutput,'all'), Parameters);
        end
    end
    
    %======================================================================
    % Helpers
    %======================================================================
    methods (Access = private)
        function Feature = generateFeatureFromInput(this, InputData, BatchSize)
            % generate features (with size D x B) by evaluating basis 
            % function on input data
            
            if BatchSize > 1
                % not support sequence dims, batch data always have batch
                % dim at ndims
                % NOTE: rely on InputData due to inconsistent output size 
                % of different model (dlarray)
                BatchDim = cellfun(@ndims, InputData,'UniformOutput',false);
                Feature = dlarray(zeros(this.FeatureSizes_(1),BatchSize),'CB');
                for ct = 1:BatchSize
                    SingleBatchData = rl.internal.dataTransformation.generalSubref(InputData,ct,BatchDim);
                    Feature(:,ct) = feval(this.BasisFcn_,SingleBatchData{:});
                end
            else
                Feature = dlarray(feval(this.BasisFcn_,InputData{:}),'CB');
            end
        end
        
        function Val = forward(this,Feature,Parameters)
            % forward pass from feature to model output. Apply post
            % processing routine determined by PostProcessFcn_
            
            % NOTE: transpose parameters because fullyconnect require size(W) = [nOut, nB]
            Parameters = Parameters';
            Val = fullyconnect(Feature,Parameters,this.Bias_);
            
            if ~isempty(this.PostProcessFcn_)
                Val = feval(this.PostProcessFcn_,Val);
            end
        end
    end
end

%======================================================================
% Helper local functions
%======================================================================
function Data = iCastDLArray(Data)
Data = cellfun(@dlarray,Data,'UniformOutput',false);
end