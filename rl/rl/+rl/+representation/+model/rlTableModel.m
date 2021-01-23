classdef rlTableModel < rl.representation.model.rlAbstractModel
    % rlTableModel Implements rlAbstractModel interface. rlTableModel is the
    %   wrapper of table.
    
    %   Copyright 2019 The MathWorks, Inc.
    
    properties (Access = private)
        InternalTable
    end
    
    methods
        function this = rlTableModel(Model,UseDevice)
            % TableModel = rlTableModel(Table,UseDevice)
            % Output
            %   - TableModel: rlTableModel object, subclass of rlAbstractModel
            % Inputs
            %   - Table     : rlTable object
            %   - UseDevice : 'cpu' or 'gpu'
            
            validateattributes(Model,{'rlTable'},{'scalar'})
            this.InternalTable = Model;
            this = updateExecutionEnvironmentImpl(this, UseDevice);
        end
    end
    
    %======================================================================
    % Implementation of abstract method
    %======================================================================
    methods (Access = protected)
        function [Val, State] = evaluateImpl(this, InputData, BatchSize, SequenceLength)
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
            %       representation level. BATCHSIZE and SEQUENCELENGTH are
            %       not use in table model.
            %
            %       STATE is an empty cell (table model does not have
            %       state)
            
            Val{1} = evaluate(this.InternalTable,InputData);
            State = {};
        end
        
        function Grad = gradientImpl(this, GradType, InputData, BatchSize, SequenceLength, LossVariable)
            % GRADIENTIMPL compute gradients from the model and input
            % data.
            %   GRAD = GRADIENTIMPL(MODEL, GRADTYPE, DATA, BATCHSIZE, SEQUENCELENGTH)
            %   computes gradients from model MODEL and input data DATA. The
            %   gradient type is determined by GradType.
            %       GRADTYPE can be 'loss-parameters', 'output-input' or 'output-parameters'
            %
            %   GRAD = GRADIENTIMPL(MODEL, 'loss-parameters', DATA, BATCHSIZE, SEQUENCELENGTH, LOSSVARIABLE)
            %   computes the gradients of the model's loss with
            %   respect to its learnable parameters. Table model can only
            %   accept MSE loss.
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
            
            % REVISIT: allow initial gradient for output-input?
            
            switch GradType
                case 'loss-parameters'
                    Grad{1} = gradientOfLossWithRespectToParameters(this.InternalTable,InputData,LossVariable);
                case 'output-input'
                    error(message('rl:agent:errUnsupportedGradient','rlTable'))
                case 'output-parameters'
                    Grad{1} = gradientWithRespectToParameters(this.InternalTable,InputData);
                otherwise
                    error(message('rl:agent:errInvalidSyntaxGradient'))
            end
        end
        
        function this = setLossImpl(this, Loss)
            % Loss is scalar function handle. Table model can only accept
            % mse loss
            
            if ~isequal(Loss, rl.util.getDefaultValueRepLoss())
                error(message('rl:agent:errTableLossCannotSet'))
            end
        end
        
        function LearnableParameters = getLearnableParametersImpl(this)
            % Return the learnable parameters of the table in a cell array
            
            LearnableParameters = {this.InternalTable.Table};
        end
        
        function this = setLearnableParametersImpl(this, LearnableParameters)
            % Set the learnable parameters of the network from a cell array
            
            this.InternalTable.Table = LearnableParameters{1};
        end
        
        function State = getStateImpl(~)
            % Return the current state of the network
            % Return empty cell if network is not RNN
            State = {};
        end
        
        function this = setStateImpl(this, ~)
            % Set the current state of the model
            % No Op
        end
        
        function this = resetStateImpl(this)
            % Set the current state of the underlying model to zeros.
            % No Op
        end
        
        function sz = getSizeImpl(this,identifier)
            % GETSIZEIMPL get the model I/O or parameter size for validation.
            % Identifier can be 'input', 'output', 'parameter'
        
            switch identifier
                case 'output'
                    sz{1} = [1,1]; % size(this.InternalTable.Table,2);
                case 'input'
                    [oinfo,ainfo] = getInfo(this.InternalTable);
                    if isempty(ainfo)
                        % V(o)
                        sz = {oinfo.Dimension};
                    else
                        % Q(o,a)
                        sz = {oinfo.Dimension ainfo.Dimension};
                    end
                case 'parameter'
                    sz = cellfun(@(x) size(x),getLearnableParameters(this),'UniformOutput',false);
            end
        end
        
        function argStruct = generateEvaluateFunctionImpl(this,argStruct)
            % Generate the evaluate function for the table model
            outputstr = argStruct.EvaluateOutputString;
            inputstr  = argStruct.EvaluateInputString;
            policyname = argStruct.PolicyName;
            matfilename = argStruct.MATFileName;
            localfcnstr = argStruct.LocalFunctionString;
            evalfcnname = argStruct.EvaluateFunctionName;
            outputstr_y = rl.codegen.handleMultipleOutputStrings(outputstr);
            
            % extract the table
            tbl = this.InternalTable.Table; %#ok<NASGU>
            
            % save the table to disk
            eval(sprintf('%s = tbl;',policyname));
            save(matfilename,policyname);
            
            % generate the local function
            load2persistentstr = sprintf([...
                'persistent %s\n',...
                'if isempty(%s)\n',...
                '\ts = coder.load(''%s'',''%s'');\n',...
                '\t%s = s.%s;\n',...
                'end'],...
                policyname,policyname,matfilename,policyname,policyname,policyname);
            [oinfo,ainfo] = getInfo(this.InternalTable);
            [oelementstr,~,~] = rl.codegen.generateFiniteElementStrings(oinfo);
            [aelementstr,~,~] = rl.codegen.generateFiniteElementStrings(ainfo);
            localevalfcnstr = sprintf([...
                'function %s = %s(%s)\n',...
                '%s\n',...
                'actionSet = %s;\n',...
                'observationSet = %s;\n',...
                'actionIndex = rl.codegen.getElementIndex(actionSet,action);\n',...
                'observationIndex = rl.codegen.getElementIndex(observationSet,%s);\n',...
                '%s = %s(observationIndex,actionIndex);\n',...
                'end'],...
                outputstr_y,evalfcnname,inputstr,...
                load2persistentstr,...
                aelementstr,...
                oelementstr,...
                argStruct.InputString,...
                outputstr_y,policyname);
            
            % attach it to the arg struct
            if isempty(localfcnstr)
                localfcnstr = localevalfcnstr;
            else
                localfcnstr = sprintf('%s\n%s',localfcnstr,localevalfcnstr);
            end
            argStruct.LocalFunctionString = localfcnstr;
        end
        
        function this = updateExecutionEnvironmentImpl(this, UseDevice)
            % UPDATEEXCUTIONENVIRONMENTIMPL updates model to use CPU or GPU
            % based on UseDevice
            
            if strcmpi(UseDevice, 'gpu')
                this.InternalTable.Table = gpuArray(this.InternalTable.Table);
            else
                if isa(this.InternalTable.Table,'gpuArray')
                    this.InternalTable.Table = gather(this.InternalTable.Table);
                end
            end
        end
        
        function HasState = hasStateImpl(this) %#ok<MANU>
            % table representation does not have state
            
            HasState = false;
        end
        
        function InternalModel = getInternalModelImpl(this)
            % GETINTERNALMODELIMPL returns the internal model
            
            InternalModel = this.InternalTable;
        end
        
        function this = setOutputTypeImpl(this,OutputType,varargin) %#ok<INUSD>
            % MODEL = SETOUTPUTTYPEIMPL(MODEL, OUTPUTTYPE) set model output to
            % be compatible with required output type
            %   OUTPUTTYPE: 'probability' or 'gaussian'
            
            error(message('rl:agent:errTableStochasticOutputNotSupport'))
        end
        
        function ExternalModel = getExternalModelImpl(this)
            % GETEXTERNALMODELIMPL returns rlTable
            
            ExternalModel = getInternalModel(this);
        end
    end
end

