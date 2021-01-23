%rlLayerModel Implements rlAbstractModel interface. rlLayerModel is the
%wrapper of neural network created from layer API.

%   Copyright 2019 The MathWorks, Inc.

classdef rlLayerModel < rl.representation.model.rlAbstractModel
    
    properties (Access = private)
        % DLT internal utilities
        InternalNetwork
        ExecutionSettings
        Assembler
        AnalyzedLayers
        Precision
        NetworkInfo
        StatefulIdx
        UseDevice
        IsStateful
        
        % Cache network I/O size
        NetworkInputSize
        
        % Nameing order for channels of the representation
        ObservationNames
        ActionNames
        
        % Index mapping from ObservationInfo and ObservationNames
        % [Obs1, Obs2, Obs3] is the observation info. This is the same order
        % as the exp comes in, then if ObservationNames was
        % [Obs2,Obs3,Obs1] then Observation index would be [2,3,1]
        ObservationIndex
        % Index mapping from ActionInfo and ActionNames
        ActionIndex
        
        % Reference for input order mapping
        InputIndex
        
        InputNames
        OutputNames
    end
    
    methods
        function this = rlLayerModel(Net, UseDevice, ObservationNames, ActionNames)
            % LayerModel = rlLayerModel(Net,UseDevice,ObservationNames,ActionNames)
            % Output
            %   - LayerModel: rlLayerModel object, subclass of rlAbstractModel
            % Inputs
            %   - Net              : neural network (layerGraph object)
            %   - UseDevoce        : 'cpu' or 'gpu'
            %   - ObservationNames : 1xN cell array of char array
            %   - ActionNames      : 1xN cell array of char array
            
            this.InternalNetwork = Net;
            this.UseDevice = UseDevice;
            this.ObservationNames = ObservationNames(:)';
            this.ActionNames = ActionNames(:)';
            
            % initial checking and replace any loss coming from the input
            % with rl FcnLossLayer (place holder)
            [this,newPlaceHolderLossLayer,lossInputName] = preBuildCheck(this);
            this = addPlaceHolderLoss(this,newPlaceHolderLossLayer,lossInputName);
            
            % build DLT internal network
            this = buildNetwork(this);
            
            % cache data size from network input/output layers
            this = cacheNetworkSize(this);
            
            % check numel(obsInfo) == numel(obsName), same for action
            % validate if all observations are inputs to representation
            this = validateObservationInput(this);
            
            % check if actions are inputs/output/none to the model
            this = checkActionIOType(this);
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
            %       representation level.
            
            % input data mapping
            InputData = this.mapInputPort(InputData);
            % cast to gpuArray if requested
            InputData = this.processInputData(InputData);
            % reshape input to correct dimension format before predict
            InputData = this.reshapeInputDimension(InputData, BatchSize, SequenceLength);
            
            if hasState(this)
                % 3rd argument is propagateState flag. Set 'true' to output
                % hidden state. In the next prediction, this hidden state
                % is passed using setState().
                % When evaluate batch of sequences input, each batch starts
                % from initial hidden state (all zeros).
                [Val, State] = statefulPredict(this.InternalNetwork,InputData,true);
            else
                Val = predict(this.InternalNetwork,InputData);
                Val = reshape(Val, 1, []);
                State = cell(size(this.InternalNetwork.Layers));
            end
            
            % REVISIT: should we unwrap gpuArray after evaluate, will make
            % compute gradient slower if the loss function calls evaluate
            Val = this.processOutputData(Val);
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
            
            % REVISIT: allow initial gradient for output-input?
            
            % input data mapping
            InputData = this.mapInputPort(InputData);
            % cast to gpuArray if requested
            % REVISIT: when to cast LossVariable to gpuArray?
            InputData = this.processInputData(InputData);
            % reshape input to correct dimension format before predict
            InputData = this.reshapeInputDimension(InputData, BatchSize, SequenceLength);
            
            switch GradType
                case 'loss-parameters'
                    Grad = this.InternalNetwork.computeGradientsForTraining(InputData,LossVariable,false);
                case 'output-input'
                    for ct = 1:numel(this.InternalNetwork.OutputSizes)
                        DummyTarget{ct} = ones([this.InternalNetwork.OutputSizes{ct} BatchSize SequenceLength],'single'); %#ok<AGROW>
                    end
                    [~, ~, ~, Grad] = this.InternalNetwork.computeGradientsForTraining( ...
                        InputData, [], false, DummyTarget);
                case 'output-parameters'
                    % LossVariable is initial gradient for this case
                    if nargin < 6
                        % if initial gradient is not provided, create
                        % adjoints of all ones
                        for ct = 1:numel(this.InternalNetwork.OutputSizes)
                            LossVariable{ct} = ones([this.InternalNetwork.OutputSizes{ct} BatchSize SequenceLength],'single');
                        end
                    end
                    Grad = this.InternalNetwork.computeGradientsForTraining(InputData,[],false,LossVariable);
                otherwise
                    error(message('rl:agent:errLayerGradientType'));
            end
        end
        
        function this = setLossImpl(this, Loss)
            % MODEL = SETLOSSIMPL(MODEL, LOSS) sets loss function handle LOSS to all
            % channel of the model MODEL. LOSS is a scalar function
            % handle.
            
            % REVISIT: support loss mapping
            % REVISIT: support multi-input loss
            
            if this.hasExternalLossLayers
                % this should not happen since any external loss are
                % replaced by place holder loss at construction of the
                % layer model
                error(message('rl:agent:errNotSupportCustomLossLayer'))
            else
                Layers  = this.InternalNetwork.OutputLayers;
                for ct = 1:numel(Layers)
                    ExternalLayer = Layers{ct}.ExternalCustomLayer;
                    ExternalLayer.LossFcn = Loss;
                    ExternalLayer.IsNetworkStateful = hasState(this);
                    this.InternalNetwork.OutputLayers{ct} = nnet.internal.cnn.layer.CustomRegressionLayer( ...
                        ExternalLayer, this.InternalNetwork.OutputLayers{ct}.LayerVerifier);
                    % REVISIT: The AnalyzedLayers is not updated with the
                    % new info from loss. Currently codegen does not
                    % support training so no issue.
                end
            end
        end
        
        function LearnableParameters = getLearnableParametersImpl(this)
            % GETLEARNABLEPARAMETERSIMPL return the learnable parameters of the
            % model.
            %
            %   LEARNABLEPARAMETERS = GETLEARNABLEPARAMETERSIMPL(MODEL) returns
            %   the learnable parameters LEARNABLEPARAMETERS of the
            %   model MODEL.
            %       LEARNABLEPARAMETERS must satisfy SETLEARNABLEPARAMETERSIMPL
            
            DLTInternalParams = this.InternalNetwork.LearnableParameters;
            LearnableParameters = cell(1,numel(DLTInternalParams));
            for ct=1:numel(DLTInternalParams)
                LearnableParameters{1,ct} = DLTInternalParams(ct).Value;
            end
        end
        
        function this = setLearnableParametersImpl(this, LearnableParameters)
            % SETLEARNABLEPARAMETERSIMPL sets the learnable parameters of the
            % model.
            %
            %   MODEL = SETLEARNABLEPARAMETERSIMPL(MODEL, LEARNABLEPARAMETERS) sets
            %   the learnable parameters of the model MODEL to
            %   LEARNABLEPARAMETERS value.
            
            % REVISIT: enable validation since setLearnableParameters
            % validateLearnableParameters(this, LearnableParameters)
            this.InternalNetwork = setLearnableParameterValues(this.InternalNetwork,LearnableParameters);
        end
        
        function State = getStateImpl(this)
            % GETSTATEIMPL return the current state of the model. If
            % the model is not stateful, return empty.
            %
            %   STATE = GETSTATEIMPL(MODEL) returns the current state STATE of
            %   the model MODEL. E.g. current hidden state of a
            %   recurrent neural network model.
            %       Special case for non stateful layer model, state is an
            %       1xNumLayer empty cell array.
            
            % Return empty cells if network is not RNN
            % REVISIT: This method currently only supports LSTM.
            % Request DLT for get state utility on internal network API
            
            State = cell(size(this.InternalNetwork.Layers));
            if hasState(this)
                for currentLayer = this.StatefulIdx
                    State{currentLayer} = ...
                        {this.InternalNetwork.Layers{currentLayer}.CellState.Value; ...
                        this.InternalNetwork.Layers{currentLayer}.HiddenState.Value};
                end
            end
        end
        
        function this = setStateImpl(this, State)
            % Set the current state of the network
            
            this.InternalNetwork = updateNetworkState(this.InternalNetwork, State);
        end
        
        function this = resetStateImpl(this)
            % Set the current state of the underlying model to zeros.
            
            this.InternalNetwork = resetNetworkState(this.InternalNetwork);
        end
        
        function sz = getSizeImpl(this,identifier)
            % GETSIZEIMPL get the model I/O or parameter size for validation.
            % Identifier can be 'input', 'output', 'parameter'
            
            switch identifier
                case 'output'
                    sz = this.InternalNetwork.OutputSizes;
                    % TODO: RNN
                    for ct = 1:numel(sz)
                        outLayerSz = sz{ct};
                        if (numel(outLayerSz) == 3 && isequal(outLayerSz(1:2), [1 1])) || ...
                                isequal(outLayerSz,1)
                            % fully connected output
                            sz{1,ct} = [outLayerSz(end) 1];
                        else
                            sz{1,ct} = outLayerSz;
                        end
                    end
                case 'input'
                    sz = this.InternalNetwork.InputSizes;
                    sz = sz(this.InputIndex);
                case 'parameter'
                    sz = cellfun(@(x) size(x),getLearnableParametersImpl(this),'UniformOutput',false);
            end
        end
        
        function argStruct = generateEvaluateFunctionImpl(this,argStruct)
            % Generate the evaluate function for the layer model
            
            outputstr = argStruct.EvaluateOutputString;
            inputstr  = argStruct.EvaluateInputString;
            policyname = argStruct.PolicyName;
            matfilename = argStruct.MATFileName;
            localfcnstr = argStruct.LocalFunctionString;
            evalfcnname = argStruct.EvaluateFunctionName;
            
            % only support SIMO network code generation
            if ~isscalar(this.InternalNetwork.InputSizes)
                error(message('rl:general:RepLayerCodeGenMIMONotSupported'));
            end
            
            % extract the DAG
            dag = iGetExternalNetwork(this.InternalNetwork,this.Assembler,this.AnalyzedLayers); %#ok<NASGU>
            
            % save the network to disk
            eval(sprintf('%s = dag;',policyname));
            save(matfilename,policyname);
            
            % generate the local evaluate network function
            load2persistentstr = sprintf([...
                'persistent %s\n',...
                'if isempty(%s)\n',...
                '\t%s = coder.loadDeepLearningNetwork(''%s'',''%s'');\n',...
                'end'],...
                policyname,policyname,policyname,matfilename,policyname);
            localevalfcnstr = sprintf([...
                'function %s = %s(%s)\n',...
                '%s\n',...
                '%s = predict(%s,%s);\n',...
                'end'],...
                rl.codegen.handleMultipleOutputStrings(outputstr),...
                evalfcnname,...
                inputstr,load2persistentstr,...
                rl.codegen.handleMultipleOutputStrings(outputstr),...
                policyname,inputstr);
            
            % attach it to the arg struct
            if isempty(localfcnstr)
                localfcnstr = localevalfcnstr;
            else
                localfcnstr = sprintf('%s\n%s',localfcnstr,localevalfcnstr);
            end
            argStruct.LocalFunctionString = localfcnstr;
            
        end
        
        function this = updateExecutionEnvironmentImpl(this,UseDevice)
            % Convert learnable params to training mode with correct
            % excecution environment (e.g. convert to gpuArray for GPU usage)
            
            this.ExecutionSettings.executionEnvironment = UseDevice;
            this.ExecutionSettings.useParallel = false;
            this.InternalNetwork = this.InternalNetwork.prepareNetworkForTraining(this.ExecutionSettings);
            this.InternalNetwork = this.InternalNetwork.optimizeNetworkForTraining(nnet.internal.cnn.optimizer.DefaultNetworkOptimizer());
            this.UseDevice = UseDevice;
        end
        
        function HasState = hasStateImpl(this)
            % HASSTATEIMPL returns true if the model has state (e.g. RNN model)
            HasState = this.IsStateful;
        end
        
        function InternalModel = getInternalModelImpl(this)
            % GETINTERNALMODELIMPL returns the internal DAGNetwork
            InternalModel = this.InternalNetwork;
        end
        
        function this = setOutputTypeImpl(this,OutputType,varargin)
            % MODEL = SETOUTPUTTYPEIMPL(MODEL, OUTPUTTYPE) set model output to
            % be compatible with required output type
            %   OUTPUTTYPE: 'probability' or 'gaussian'
            % e.g. MODEL = SETOUTPUTTYPEIMPL(MODEL, 'probability') add
            % softmax to the model output to ensure it outputs probability
            
            % REVISIT: only support adding softmax for probability output
            % type. For now, 'gaussian' output must be embedded in the
            % neural network (Layer API)
            if strcmpi(OutputType,'probability')
                % REVISIT: we dont need to convert the model back to
                % layerGraph to check for softmax layer. This will
                % improve performance if input model already has
                % softmax layer. We can use OutputNames property.
                
                DAG = iGetExternalNetwork(this.InternalNetwork,this.Assembler,this.AnalyzedLayers);
                LGraph = layerGraph(DAG);
                [DiGraph,LayerNodes] = iGetDigraph(LGraph);
                % then get the adjacency matrix
                AdjMatrix = DiGraph.adjacency';
                % find the output (sink) nodes by finding nodes that are
                % NOT inputs to any other node.
                LossLayerIdx = find(~any(AdjMatrix,1));
                % get the layer names
                Names = {LGraph.Layers.Name};
                
                DoAttachSoftmax = false;
                for ct = 1:numel(LossLayerIdx)
                    yi = LossLayerIdx(ct);
                    LossLayer = LayerNodes(yi);
                    % output layer is the predecessor of the loss layer
                    OutputLayer = LayerNodes(predecessors(DiGraph,yi));
                    if ~ (isa(OutputLayer,'nnet.internal.cnn.layer.SoftmaxLayer') || isa(OutputLayer,'nnet.cnn.layer.SoftmaxLayer'))
                        DoAttachSoftmax = true;
                        Names = matlab.lang.makeUniqueStrings([Names {'RepresentationSoftMax'}]);
                        OutputLayerNotSoftmax{ct} = OutputLayer.Name; %#ok<AGROW>
                        LossLayersToRemove{ct} = LossLayer.Name; %#ok<AGROW>
                        LossLayersToAttach{ct} = [softmaxLayer('Name',Names{end}) LossLayer]; %#ok<AGROW>
                        % REVISIT: update output name to softmax layer name
                        % not necessary now since stochastic actor not
                        % require naming
                    end
                end
                
                if DoAttachSoftmax
                    % remove all loss layer from the layer graph
                    LGraph = removeLayers(LGraph,LossLayersToRemove);
                    for ct = 1:numel(OutputLayerNotSoftmax)
                        % attach softmax layers + loss layers
                        LGraph = addLayers(LGraph,LossLayersToAttach{ct});
                        % connect networks outputs to input of softmax
                        LossInputName = LossLayersToAttach{ct}(1).Name;
                        LGraph = connectLayers(LGraph,OutputLayerNotSoftmax{ct},LossInputName);
                    end
                    this.InternalNetwork = LGraph;
                    % build DLT internal network
                    this = buildNetwork(this);
                end
            end
        end
        
        function ExternalModel = getExternalModelImpl(this)
            % EXTERNALMODEL = GETEXTERNALMODEL(MODEL) returns the external
            % SeriesNetwork or DAGNetwork
            % The RL Custom loss layer is replace with regressionLayer or
            % classificationLayer (only if network has softmaxLayer at output)
            
            % create layerGraph from internal network
            LGraph = layerGraph(iGetExternalNetwork(this.InternalNetwork,this.Assembler,this.AnalyzedLayers));
            
            % replace RL loss layer with dummy loss layer
            % NOTE: only support single output network for now
            if numel(LGraph.OutputNames) > 1
                error(message('rl:agent:errLayerMultiChannelGetExternalNotSupport'))
            end
            LossLayerName = LGraph.OutputNames{1};
            OutputLayerIdx = arrayfun(@(x) strcmpi(x.Name,this.OutputNames{1}),LGraph.Layers);
            OutputLayer = LGraph.Layers(OutputLayerIdx);
            UseClassificationLayer = isa(OutputLayer,'nnet.cnn.layer.SoftmaxLayer');
            if UseClassificationLayer
                DummyLayer = classificationLayer('Name',LossLayerName,'Classes','auto');
            else
                DummyLayer = regressionLayer('Name',LossLayerName);
            end
            LGraph = replaceLayer(LGraph,LossLayerName,DummyLayer);
            
            % construct external DAG network from layerGraph
            if UseClassificationLayer
                % temporary bypass DLT warning on classification layer does
                % not have classes property
                warning('off','nnet_cnn:internal:cnn:analyzer:NetworkAnalyzer:NetworkHasWarnings')
                ExternalModel = assembleNetwork(LGraph);
                warning('on','nnet_cnn:internal:cnn:analyzer:NetworkAnalyzer:NetworkHasWarnings')
            else
                ExternalModel = assembleNetwork(LGraph);
            end
        end
    end
    
    %======================================================================
    % DLT Internal Utilities
    %======================================================================
    methods (Access = private)
        function [this,newPlaceHolderLossLayer,lossInputName] = preBuildCheck(this)
            % Check network validity and cache some utility properties
            % If the input layerGraph has any loss layer, return new
            % place holder layer
            % (with loss translated from external layer and same name with the old layer)
            % and name of the input layer to the loss
            % NOTE:
            % - This assumes loss layer can only have single inputs
            % - Throw an error if detect custom output layer. Suggest use
            % setLoss instead.
            
            [dg,layerNodes] = iGetDigraph(this.InternalNetwork);
            
            % check for any unsupported layers
            unsupportedLayers = iGetUnsupportedLayer();
            for ct = 1:numel(unsupportedLayers)
                unsupportedSingleLayer = unsupportedLayers(ct);
                if any(arrayfun(@(x) isa(x,unsupportedSingleLayer),layerNodes))
                    error(message('rl:agent:errInvalidLayers'));
                end
            end
            
            % TODO: why do we make sure that Normalization is set to None in imageInputLayer
            inputLayer = layerNodes(arrayfun(@(x) isa(x,'nnet.cnn.layer.ImageInputLayer'),layerNodes));
            if any(arrayfun(@(x) ~strcmpi(x.Normalization,'none'),inputLayer))
                error(message('rl:agent:errImageInputLayerNormalizationNone'));
            end
            
            % make sure that OutputMode is set to sequence in lstmLayer
            inputLayer = layerNodes(arrayfun(@(x) isa(x,'nnet.cnn.layer.LSTMLayer'),layerNodes));
            if any(arrayfun(@(x) ~strcmpi(x.OutputMode,'sequence'),inputLayer))
                error(message('rl:agent:errLSTMLayerOutputSequence'));
            end
            
            % make sure the IOs are defined
            % get the adjacency matrix
            adj = dg.adjacency';
            % find the output (sink) nodes by finding nodes that are
            % NOT inputs to any other node.
            yidx = find(~any(adj,1));
            
            % detect and remove any loss layer from layerGraph
            % attempt to translate external DLT loss layer into rl custom
            % loss
            newPlaceHolderLossLayer = [];
            lossInputName = {};
            for i = 1:numel(yidx)
                yi = yidx(i);
                outLayer = layerNodes(yi);
                if iIsaLossLayer(outLayer)
                    % if model includes a loss layer, the output layer is the
                    % predecessor of the loss
                    yidx(i) = predecessors(dg,yi);
                    
                    % if model includes a DLT layer, translate to RL loss
                    % function and replace with RL loss layer
                    outLayerIsDLTLossLayer = iIsaDLTLossLayer(outLayer);
                    if outLayerIsDLTLossLayer
                        if isa(outLayer,'nnet.cnn.layer.RegressionOutputLayer')
                            % regressionLayer
                            loss = @rl.loss.rlmse;
                        elseif isa(outLayer,'nnet.cnn.layer.ClassificationOutputLayer') || isa(outLayer,'rl.layer.internal.ClassificationLayerWithEntropy')
                            % classificationLayer
                            loss = @rl.loss.policyGradientDiscrete;
                        elseif isa(outLayer,'rl.layer.internal.ClipPPOLossLayer')
                            % previous releases only have discrete PPO
                            loss = @rl.loss.ppoClipDiscrete;
                        else
                            error(message('rl:agent:errNotSupportCustomLossLayer'));
                        end
                        % track the name of the output layer to be connect to
                        % FcnLossLayer
                        lossInputName = cat(1,lossInputName,{layerNodes(yidx(i)).Name});
                        % create the FcnLossLayer with the loss translated from
                        % external DLT loss layer
                        newPlaceHolderLossLayer = cat(1,newPlaceHolderLossLayer,...
                            rl.layer.FcnLossLayer('Name',outLayer.Name,'LossFcn',loss));
                    end
                end
            end
            
            % find input (src) nodes by finding nodes that have NO incoming
            % connections
            uidx = ~any(adj,2)';
            
            % get the io names
            fcn = @(x) x.Name;
            yNames = arrayfun(fcn,layerNodes(yidx),'UniformOutput',false);
            uNames = arrayfun(fcn,layerNodes(uidx),'UniformOutput',false);
            
            % cache input nodes, output nodes
            this.InputNames  = uNames(:)';
            this.OutputNames = yNames(:)';
        end
        
        function this = addPlaceHolderLoss(this,newPlaceHolderLossLayer,lossInputName)
            % Add rl FcnLossLayer to the layer graph so an internal network
            % can be built.
            %   - If newPlaceHolderLossLayer is empty, attach loss with
            %   empty forward method (place holder).
            %   - If newPlaceHolderLossLayer is not empty, remove existing
            %   loss layers from the layerGraph and attach
            %   newPlaceHolderLossLayer to lossInputName (output layers)
            %   newPlaceHolderLossLayer is FcnLossLayer with predefined loss
            % NOTE: this assumes loss layer can only have single inputs
            
            if isempty(newPlaceHolderLossLayer)
                % get the digraph to the network
                dg = iGetDigraph(this.InternalNetwork);
                % then get the adjacency matrix
                adj = dg.adjacency';
                % find the output (sink) nodes by finding nodes that are
                % NOT inputs to any other node.
                snkidx = ~any(adj,1);
                % get the layer names
                Names = {this.InternalNetwork.Layers.Name};
                % get the output layer name
                OutputLayerName = Names(snkidx);
                OutputLayer = this.InternalNetwork.Layers(snkidx);
                
                % add place holder loss to each output layer only if
                % internal loss layer is not already attached
                InternalLossLayerIdx = arrayfun(@(x) isa(x,'rl.layer.FcnLossLayer'), OutputLayer);
                InternalLossLayerName = OutputLayerName(InternalLossLayerIdx);
                OutputLayerName(InternalLossLayerIdx) = [];
                NumOutputLayer = numel(OutputLayerName);
                
                if NumOutputLayer > 0
                    LossName = matlab.lang.makeUniqueStrings([repmat({'RepresentationLoss'},1,NumOutputLayer) InternalLossLayerName]);
                    tempRepresentation = this.InternalNetwork;
                    for ct = 1:NumOutputLayer
                        LossLayer = rl.layer.FcnLossLayer(...
                            'Name',LossName{ct});
                        tempRepresentation = addLayers(tempRepresentation,LossLayer);
                        tempRepresentation = connectLayers(tempRepresentation,OutputLayerName{ct},LossName{ct});
                    end
                    
                    this.InternalNetwork = tempRepresentation;
                end
            else
                % remove existing loss (same name with newPlaceHolderLossLayer)
                this.InternalNetwork = removeLayers(this.InternalNetwork, {newPlaceHolderLossLayer.Name});
                % connect new FcnLossLayer
                this.InternalNetwork = addLayers(this.InternalNetwork, newPlaceHolderLossLayer);
                for ct = numel(newPlaceHolderLossLayer)
                    this.InternalNetwork = connectLayers(this.InternalNetwork, lossInputName{ct}, newPlaceHolderLossLayer(ct).Name);
                end
            end
        end
        
        function this = buildNetwork(this)
            % Build internal DLT network only when loss is attached
            
            % Set desired precision
            this.Precision = nnet.internal.cnn.util.Precision('single');
            
            [this.InternalNetwork, this.ExecutionSettings, ...
                this.Assembler, this.AnalyzedLayers,this.NetworkInfo] = createInternalNeuralNetwork(this);
            
            % Stateful utility
            this.IsStateful = any(nnet.internal.cnn.util.isStatefulLayer( this.InternalNetwork.Layers ));
            this.StatefulIdx = find( nnet.internal.cnn.util.isStatefulLayer(this.InternalNetwork.Layers) );
        end
        
        function [internalNetwork,executionSettings,assembler,analyzedLayers,networkInfo] = createInternalNeuralNetwork(this)
            lgraph = this.InternalNetwork;
            
            % Check network for SISO & MIMO
            isaDAG = isa(lgraph,'nnet.cnn.LayerGraph');
            analyzedLayers = nnet.internal.cnn.analyzer.NetworkAnalyzer(lgraph);
            analyzedLayers.applyConstraints([
                rl.util.Architecture()
                nnet.internal.cnn.analyzer.constraints.Connections()
                nnet.internal.cnn.analyzer.constraints.ConnectedComponents()
                nnet.internal.cnn.analyzer.constraints.CustomLayers()
                nnet.internal.cnn.analyzer.constraints.LSTM()
                nnet.internal.cnn.analyzer.constraints.Names()
                nnet.internal.cnn.analyzer.constraints.Propagation()]);
            analyzedLayers.throwIssuesIfAny()
            internalLayers = analyzedLayers.InternalLayers;
            networkInfo = nnet.internal.cnn.util.ComputeNetworkInfo(isaDAG,internalLayers);
            
            % Set up and validate training
            dummyDLTOption = iCreateDummyDLTOptions(this.UseDevice);
            executionSettings = nnet.internal.cnn.assembler.setupExecutionEnvironment(...
                dummyDLTOption, networkInfo.IsRNN, [], this.Precision );
            
            % Assemble internal network
            strategy = nnet.internal.cnn.assembler.NetworkAssemblerStrategyFactory.createStrategy(~networkInfo.IsDAG); % try with removing it
            assembler = nnet.internal.cnn.assembler.TrainingNetworkAssembler(strategy);
            
            internalNetwork = assembler.assemble(analyzedLayers, executionSettings);
            internalNetwork.optimizeNetworkForTraining(nnet.internal.cnn.optimizer.DefaultNetworkOptimizer());
        end
        
        function flag = hasExternalLossLayers(this)
            ix = arrayfun(@(x) iIsaDLTLossLayer(x),this.InternalNetwork.Layers);
            flag = any(ix);
        end
        
        function validateLearnableParameters(this, LearnableParameters)
            
            if ~isa(LearnableParameters,'cell')
                % TODO: is cell specific to neural network? If yes need to
                % enhance error message to be more descriptive
                error(message('rl:agent:errLearnableParameters'));
            end
            ParamSize = getSizeImpl(this,'parameter');
            ValSize = cellfun(@(x) size(x),LearnableParameters,'UniformOutput',false);
            if length(ParamSize)~=length(ValSize)
                error(message('rl:agent:errParameterSize'));
            end
            if ~all(cellfun(@(x,y) isequal(x,y),ParamSize,ValSize))
                error(message('rl:agent:errParameterSize'));
            end
        end
    end
    
    %======================================================================
    % Data manipulation
    %======================================================================
    methods (Access = private)
        function Data = mapInputPort(this,Data)
            % Map input data to correct input port order
            
            [~,ix] = sort(this.InputIndex);
            if numel(ix) ~= numel(Data)
                % handle when number of input data cells different from
                % number of specifed input channel
                error(message('rl:agent:errRepIncorrectNumInputChannel'));
            end
            Data = Data(ix);
        end
        
        function Data = reshapeInputDimension(this, Data, BatchSize, SequenceLength)
            % Reshape input data to be compatible with neural networks
            % operation
            % e.g. [d d nBatch] to [d d 1 nBatch]
            % e.g. [d d nBatch nSequence] to [d d 1 nBatch nSequence]
            
            % do not support BatchSize of data with SequenceLength = 1 for
            % RNN
            if (BatchSize > 1) && (SequenceLength == 1) && hasState(this)
                error(message('rl:agent:errLayerModelNotSupportBatchDataSequenceLengthEq1'))
            end
            for ct = 1:numel(Data)
                Data{ct} = reshape(single(Data{ct}), [this.NetworkInputSize{ct} BatchSize SequenceLength]);
            end
        end
        
        function this = cacheNetworkSize(this)
            this.NetworkInputSize = this.InternalNetwork.InputSizes;
        end
        
        function Data = processInputData(this,Data)
            % Cast input data to gpuArray if the model uses GPU
            % NOTE: casting with DLT precision object will causes issue
            % when calling rep.evaluate on the loss function since it
            % calls extractdata() if input is dlarray
            
            if isGPU(this)
                % Data = cellfun(@(x) gpuArray(this.Precision.cast(x)),Data,'UniformOutput',false);
                Data = cellfun(@(x) gpuArray(single(x)),Data,'UniformOutput',false);
            else
                % Data = cellfun(@(x) this.Precision.cast(gather(x)),Data,'UniformOutput',false);
                Data = cellfun(@(x) single(gather(x)),Data,'UniformOutput',false);
            end
        end
        
        function Data = processOutputData(this,Data)
            % Cast output data to numeric if the model uses GPU
            % REVISIT: map output channel?
            if isGPU(this)
                Data = cellfun(@(x) gather(x), Data, 'UniformOutput',false);
            end
        end
        
        function flag = isGPU(this)
            flag = ismember(this.ExecutionSettings.executionEnvironment, {'gpu'} );
        end
    end
    
    %======================================================================
    % Validation
    %======================================================================
    methods (Access = private)
        function this = validateObservationInput(this)
            % Check if all observations are inputs to representation and
            % set input port index
            
            [IntersectObsName,~,ObsIdx] = intersect(this.ObservationNames, this.InputNames, 'stable');
            ObsIdx = ObsIdx(:)';
            if isempty(IntersectObsName) || (numel(IntersectObsName) < numel(this.ObservationNames))
                error(message('rl:agent:errObservationNames'));
            end
            this.ObservationIndex = ObsIdx;
            this.InputIndex = ObsIdx;
        end
        
        function this = checkActionIOType(this)
            % Check if all actions are inputs to the representation and
            % set input port index
            
            [IntersectActionInputName,~,ActionIdx] = intersect(this.ActionNames, this.InputNames, 'stable');
            ActionIdx = ActionIdx(:)';
            if ~isempty(IntersectActionInputName)
                ActionIdx = reshape(ActionIdx,1,[]);
                this.ActionIndex = {"input", ActionIdx};
                this.InputIndex = [this.InputIndex this.ActionIndex{2}];
            end
            
            % Both cases cannot happen since DLT checks for unique layer
            % names. Would already error out before reaching this step
            [IntersectActionOutputName,~,ActionIdx] = intersect(this.ActionNames, this.OutputNames, 'stable');
            ActionIdx = ActionIdx(:)';
            if ~isempty(IntersectActionOutputName)
                this.ActionIndex = {"output", ActionIdx};
                % REVISIT: output index for output mapping?
            end
            
            if (max(numel(IntersectActionInputName),numel(IntersectActionOutputName)) < numel(this.ActionNames))
                error(message('rl:agent:errActionNames'));
            end
        end
    end
    
    %======================================================================
    % Save/Load
    %======================================================================
    methods
        function s = saveobj(this)
            % Turn internal DLT API to external DAGNetwork object when
            % saving. Layer model will be reconstructed in loadobj
            
            % TODO remove internal loss layers and replace with
            % regression
            
            % save properties required for reconstruction
            s.ObservationNames = this.ObservationNames;
            s.ActionNames = this.ActionNames;
            
            % turn network to user-facing layerGraph object
            % REVISIT: is this a GPU network or CPU? If GPU, might not able
            % to load on a non-GPU machine
            s.DLTModelExternal = iGetExternalNetwork(this.InternalNetwork,this.Assembler,this.AnalyzedLayers);
            s.DLTModelExternal = rl.internal.dataTransformation.networkToLayerGraph(s.DLTModelExternal);
            
            % save representation with UseDevice set to 'cpu'
            s.SaveData.UseDevice = this.UseDevice;
            s.UseDevice = 'cpu';
        end
    end
    
    methods (Static)
        function this = loadobj(s)
            if isstruct(s)
                % revert to gpuArray if UseDevice is 'gpu'
                s.UseDevice = s.SaveData.UseDevice;
                this = rl.representation.model.rlLayerModel(...
                    s.DLTModelExternal,s.UseDevice,s.ObservationNames,s.ActionNames);
            else
                this = s;
                % convert learnable params to training mode
                % UseDevice property is updated to s.SaveData.UseDevice
                this = updateExecutionEnvironment(this,s.SaveData.UseDevice);
            end
        end
    end
end

%==========================================================================
% Helper local functions
%==========================================================================
function DLTOpts = iCreateDummyDLTOptions(UseDevice)
DLTOpts = trainingOptions('adam','ExecutionEnvironment',UseDevice);
end

function flag = iIsaLossLayer(layer)
flag = isa(layer,'rl.layer.FcnLossLayer') || iIsaDLTLossLayer(layer);
end

function flag = iIsaDLTLossLayer(layer)
% Check if the layer is a DLT layer but not RL internal FcnLossLayer

if isa(layer,'rl.layer.FcnLossLayer')
    % include this path since FcnLossLayer is a nnet.layer.RegressionLayer
    flag = false;
else
    flag = ...
        isa(layer,'nnet.cnn.layer.RegressionOutputLayer') || ... built-in regressionLayer
        isa(layer,'nnet.cnn.layer.ClassificationOutputLayer') || ... built-in classificationLayer
        isa(layer,'nnet.layer.RegressionLayer') || ... custom regression layer
        isa(layer,'nnet.layer.ClassificationLayer') || ... custom classification layer
        isa(layer,'rl.layer.internal.ClassificationLayerWithEntropy') || ...
        isa(layer,'rl.layer.internal.ClipPPOLossLayer');
end
end

function unsupportedLayers = iGetUnsupportedLayer
unsupportedLayers = [
    "nnet.cnn.layer.image3dInputLayer"
    "nnet.cnn.layer.BatchNormalizationLayer"
    "nnet.cnn.layer.BiLSTMLayer"
    "nnet.cnn.layer.GRULayer"
    "nnet.cnn.layer.DicePixelClassificationLayer"
    "nnet.cnn.layer.FocalLossLayer"
    "nnet.cnn.layer.PixelClassificationLayer"
    "nnet.cnn.layer.RCNNBoxRegressionLayer"
    "nnet.cnn.layer.RegionProposalLayer"
    "nnet.cnn.layer.RPNClassificationLayer"
    "nnet.cnn.layer.YOLOv2OutputLayer"
    ];
end

function [dg,layerNodes] = iGetDigraph(lgraph)
% get the digraph to the layerGraph
dg = extractPrivateDirectedGraph(lgraph);
layerNodes = dg.Nodes{:,1};
end

function DAGNet = iGetExternalNetwork(InternalDAGNetwork,Assembler,AnalyzedLayers)
TempadModel = InternalDAGNetwork.prepareNetworkForPrediction();
TempadModel = TempadModel.setupNetworkForHostPrediction();
DAGNet = Assembler.createExternalNetwork(TempadModel,AnalyzedLayers);
end