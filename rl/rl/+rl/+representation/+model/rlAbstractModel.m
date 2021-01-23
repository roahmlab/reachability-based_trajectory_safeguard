%rlAbstractModel Defines the interface of computation model for
% reinforcement learning representation. It contains a Model object, which
% can be neural networks created from layer API or dlnetwork, rlTable,
% function handle.
%
% rlAbstractModel properties:
%    Model          - Computation engine, can be neural networks created
%                   from layer API or dlnetwork, rlTable, function handle.
%
% rlAbstractModel public methods:
%    evaluate - Evaluate the model from given input data.
%    gradient - Compute gradient with various types. E.g. dLoss/dParameter
%    setLoss  - Set the model's loss.
%    syncParameters - Update learnable parameters of one model
%                     with another
%    get/setLearnableParameters - Get or set the model's parameter
%    get/setState - Get or set the model's state
%    resetState - Reset the model's state
%    generateEvaluateFunction - generate codegen function

% Copyright 2019 The MathWorks, Inc.

classdef rlAbstractModel
    
    %======================================================================
    % Public API
    %======================================================================
    methods
        function [Val, State] = evaluate(this, InputData, BatchSize, SequenceLength)
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
            
            [Val, State] = evaluateImpl(this, InputData, BatchSize, SequenceLength);
        end
        
        function Grad = gradient(this, GradType, InputData, BatchSize, SequenceLength, varargin)
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
            %   GRAD = GRADIENT(MODEL, 'output-parameters', DATA, BATCHSIZE, SEQUENCELENGTH)
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
            %
            %       BATCHSIZE and SEQUENCELENGTH are scalar given by the
            %       representation level.
            %
            %   GRAD = GRADIENT(MODEL, 'output-parameters', DATA, BATCHSIZE, SEQUENCELENGTH, INITIALGRADIENT)
            %   computes the gradients of the model's outputs from given
            %   initial gradients with respect to its inputs.
            %
            %       INITIALGRADIENT must have the same data type, shape and
            %       size with the model output.
            
            % REVISIT: If LOSSVARIABLE is a cell array with numel == number
            % of model's outputs, each cell will be passed to the loss
            % function of each output respectively.
            
            Grad = gradientImpl(this, GradType, InputData, BatchSize, SequenceLength, varargin{:});
        end
        
        function this = syncParameters(this, Model, SmoothFactor)
            % SYNCPARAMETERS update the learnable parameters of one
            % model with another.
            %
            %   TARGETMODEL = SYNCPARAMETERS(TARGETMODEL, SOURCEMODEL, SMOOTHFACTOR)
            %   updates the learnable parameters of the TARGETMODEL
            %   model with thos of SOURCEMODEL.
            %   TARGETPARAMS = SMOOTHFACTOR*SOURCEPARAMS + (1 - SMOOTHFACTOR)*TARGETPARAMS
            
            NewParameterValues = getLearnableParameters(Model);
            CurrentParameterValues = getLearnableParameters(this);
            UpdateFun = @(x,y) x*(1-SmoothFactor) + y*SmoothFactor;
            this = setLearnableParameters(this,dlupdate(UpdateFun,CurrentParameterValues,NewParameterValues));
        end
        
        function this = setLoss(this, Loss)
            %MODEL = SETLOSS(MODEL, LOSS) sets loss function handle LOSS to all
            % channel of the model MODEL. LOSS is a scalar function
            % handle.
            
            this = setLossImpl(this, Loss);
        end
        
        function LearnableParameters = getLearnableParameters(this)
            % GETLEARNABLEPARAMETERS return the learnable parameters of the
            % model.
            %
            %   LEARNABLEPARAMETERS = GETLEARNABLEPARAMETERS(MODEL) returns
            %   the learnable parameters LEARNABLEPARAMETERS of the
            %   model MODEL.
            %       LEARNABLEPARAMETERS must satisfy SETLEARNABLEPARAMETERS
            
            LearnableParameters = getLearnableParametersImpl(this);
        end
        
        function this = setLearnableParameters(this, LearnableParameters)
            % SETLEARNABLEPARAMETERS sets the learnable parameters of the
            % model.
            %
            %   MODEL = SETLEARNABLEPARAMETERS(MODEL, LEARNABLEPARAMETERS) sets
            %   the learnable parameters of the model MODEL to
            %   LEARNABLEPARAMETERS value.
            
            this = setLearnableParametersImpl(this, LearnableParameters);
        end
        
        function State = getState(this)
            % GETSTATE return the current state of the model. If
            % the model is not stateful, return empty.
            %
            %   STATE = GETSTATE(MODEL) returns the current state STATE of
            %   the model MODEL. E.g. current hidden state of a
            %   recurrent neural network model.
            %       Special case for non stateful layer model, state is an
            %       1xNumLayer empty cell array.
            
            State = getStateImpl(this);
        end
        
        function this = setState(this, State)
            % SETSTATE set the current state of the model. No-op
            % for model that is not stateful.
            %
            %   MODEL = SETSTATE(MODEL, STATE) sets the current state of the
            %   model MODEL to STATE. STATE must have similar size
            %   with data coming from GETSTATE
            
            this = setStateImpl(this, State);
        end
        
        function this = resetState(this)
            % RESETSTATE reset the state of the model to zeros.
            % No-op for model that is not stateful.
            %
            %   MODEL = RESETSTATE(MODEL) resets the current state of the
            %   model MODEL to zeros.
            
            this = resetStateImpl(this);
        end
        
        function argStruct = generateEvaluateFunction(this,argStruct)
            % GENERATEPOLICYFUNCTION generate a stand-alone static policy
            % function. Useful for code generation
            
            argStruct = generateEvaluateFunctionImpl(this,argStruct);
        end
        
        function this = updateExecutionEnvironment(this,UseDevice)
            % UPDATEEXCUTIONENVIRONMENT updates model to use CPU or GPU
            % based on UseDevice
            
            this = updateExecutionEnvironmentImpl(this,UseDevice);
        end
        
        function InternalModel = getInternalModel(this)
            % GETINTERNALMODEL returns the internal model
            
            InternalModel = getInternalModelImpl(this);
        end
        
        function sz = getSize(this,Identifier)
            % GETSIZE get the model I/O or parameter size for validation.
            % Identifier can be 'input', 'output', 'parameter'
            
            sz = getSizeImpl(this,Identifier);
        end
        
        function HasState = hasState(this)
            % HASSTATE returns true if the model has state (e.g. RNN model)
            HasState = hasStateImpl(this);
        end
        
        function this = setOutputType(this,OutputType,varargin)
            % MODEL = SETOUTPUTTYPE(MODEL, OUTPUTTYPE) set model output to 
            % be compatible with required output type
            %   OUTPUTTYPE: 'probability' or 'gaussian'
            %   if OUTPUTTYPE is 'gaussian', varargin contains a struct of
            %   scale and bias
            % e.g. MODEL = SETOUTPUTTYPE(MODEL, 'probability') add
            % softmax to the model output to ensure it outputs probability
            
            OutputType = validatestring(OutputType,["probability","gaussian"], 'setOutputType', 'setOutputType', 2);
            this = setOutputTypeImpl(this,OutputType,varargin{:});
        end
        
        function ExternalModel = getExternalModel(this)
            % EXTERNALMODEL = GETEXTERNALMODEL(MODEL) returns the external
            % computation model
            % * Deep Neural Network Model
            %       EXTERNALMODEL = GETEXTERNALMODEL(MODEL) returns a 
            %       SeriesNetwork, DAGNetwork or dlnetwork object.
            % * Table Model
            %       EXTERNALMODEL = GETEXTERNALMODEL(MODEL) returns an 
            %       rlTable object.
            % * Basis Function Model
            %       EXTERNALMODEL = GETEXTERNALMODEL(MODEL) returns a cell
            %       contains the basis function and the learnable parameters.
            
            ExternalModel = getExternalModelImpl(this);
        end
    end
    
    %======================================================================
    % Protected methods
    %======================================================================
    methods (Abstract, Access = protected)
        [Val, State] = evaluateImpl(this, InputData, BatchSize, SequenceLength)
        % EVALUATEIMPL evaluate the model with given input data.
        %   [VALUE, STATE] = EVALUATEIMPL(MODEL, DATA, BATCHSIZE, SEQUENCELENGTH)
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
        
        Grad = gradientImpl(this, GradType, InputData, BatchSize, SequenceLength, varargin)
        % GRADIENTIMPL compute gradients from the model and input
        % data.
        %   GRAD = GRADIENTIMPL(MODEL, GRADTYPE, DATA, BATCHSIZE, SEQUENCELENGTH)
        %   computes gradients from model MODEL and input data DATA. The
        %   gradient type is determined by GradType.
        %       GRADTYPE can be 'loss-parameters', 'output-input' or 'output-parameters'
        %
        %   GRAD = GRADIENTIMPL(MODEL, 'loss-parameters', DATA, BATCHSIZE, SEQUENCELENGTH, LOSSVARIABLE)
        %   computes the gradients of the model's loss with
        %   respect to its learnable parameters. Must set loss with
        %   SETLOSS function to use this gradient type.
        %
        %   GRAD = GRADIENTIMPL(MODEL, 'input-output', DATA, BATCHSIZE, SEQUENCELENGTH)
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
        %   GRAD = GRADIENTIMPL(MODEL, 'output-parameters', DATA, BATCHSIZE, SEQUENCELENGTH)
        %   computes the gradients of the model's output with
        %   respect to its learnable parameters.
        
        this = setLossImpl(this, Loss)
        % MODEL = SETLOSSIMPL(MODEL, LOSS) sets loss function handle LOSS to all
        % channel of the model MODEL. LOSS is a scalar function
        % handle.
        
        LearnableParameters = getLearnableParametersImpl(this)
        % GETLEARNABLEPARAMETERSIMPL return the learnable parameters of the
        % model.
        %
        %   LEARNABLEPARAMETERS = GETLEARNABLEPARAMETERSIMPL(MODEL) returns
        %   the learnable parameters LEARNABLEPARAMETERS of the
        %   model MODEL.
        %       LEARNABLEPARAMETERS must satisfy SETLEARNABLEPARAMETERSIMPL
        
        this = setLearnableParametersImpl(this, LearnableParameters)
        % SETLEARNABLEPARAMETERSIMPL sets the learnable parameters of the
        % model.
        %
        %   MODEL = SETLEARNABLEPARAMETERSIMPL(MODEL, LEARNABLEPARAMETERS) sets
        %   the learnable parameters of the model MODEL to
        %   LEARNABLEPARAMETERS value.
        
        State = getStateImpl(this)
        % GETSTATEIMPL return the current state of the model. If
        % the model is not stateful, return empty.
        %
        %   STATE = GETSTATEIMPL(MODEL) returns the current state STATE of
        %   the model MODEL. E.g. current hidden state of a
        %   recurrent neural network model.
        %       Special case for non stateful layer model, state is an
        %       1xNumLayer empty cell array.
        
        this = setStateImpl(this, State)
        % SETSTATEIMPL set the current state of the model. No-op
        % for model that is not stateful.
        %
        %   MODEL = SETSTATEIMPL(MODEL, STATE) sets the current state of the
        %   model MODEL to STATE. STATE must have similar size
        %   with data coming from GETSTATE
        
        this = resetStateImpl(this)
        % RESETSTATEIMPL reset the state of the model to zeros.
        % No-op for model that is not stateful.
        %
        %   MODEL = RESETSTATEIMPL(MODEL) resets the current state of the
        %   model MODEL to zeros.
        
        this = updateExecutionEnvironmentImpl(this,UseDevice)
        % UPDATEEXCUTIONENVIRONMENTIMPL updates model to use CPU or GPU
        % based on UseDevice
        
        HasState = hasStateImpl(this)
        % HASSTATEIMPL returns true if the model has state (e.g. RNN model)
        
        sz = getSizeImpl(this,Identifier)
        % GETSIZEIMPL get the model I/O or parameter size for validation.
        % Identifier can be 'input', 'output', 'parameter'
        
        InternalModel = getInternalModelImpl(this);
        % GETINTERNALMODELIMPL returns the internal model
        
        this = setOutputTypeImpl(this,OutputType,varargin);
        % MODEL = SETOUTPUTTYPEIMPL(MODEL, OUTPUTTYPE) set model output to
        % be compatible with required output type
        %   OUTPUTTYPE: 'probability' or 'gaussian'
        % e.g. MODEL = SETOUTPUTTYPEIMPL(MODEL, 'probability') add
        % softmax to the model output to ensure it outputs probability
        
        ExternalModel = getExternalModelImpl(this)
        % EXTERNALMODEL = GETEXTERNALMODELIMPL(MODEL) returns the external
        % computation model
        % * Deep Neural Network Model
        %       EXTERNALMODEL = GETEXTERNALMODELIMPL(MODEL) returns a SeriesNetwork, DAGNetwork or
        %       dlnetwork object.
        % * Table Model
        %       EXTERNALMODEL = GETEXTERNALMODELIMPL(MODEL) returns an rlTable object.
        % * Basis Function Model
        %       EXTERNALMODEL = GETEXTERNALMODELIMPL(MODEL) returns a cell contains the basis
        %       function and the learnable parameters.
    end
    
    methods (Access = protected)
        function argStruct = generateEvaluateFunctionImpl(this,argStruct) %#ok<INUSD>
            % model must overload this method for code gen support
            error(message('rl:general:RepCodeGenNotSupported',class(this)));
        end
    end
end