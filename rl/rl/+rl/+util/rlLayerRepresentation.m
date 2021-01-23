classdef rlLayerRepresentation < rl.util.rlAbstractRepresentation
% rlLayerRepresentation: RL representation for RL agents using a neural network
% from Deep Learning Toolbox
%
% To create an rlLayerRepresention object, use the rlRepresentation
% function.
%
%   rep = rlRepresentation(ADMODEL,OINFO,'Observation',ONAMES) creates a
%   representation with default options using ADMODEL, which is a neural
%   network from the Deep Learning Toolbox. Specify the network input layer names,
%   ONAMES, associated with each observation specification as a cell array.
%   The names in ONAMES must be the same order as the observation
%   specifications, OINFO.
%
%   rep = rlRepresentation(ADMODEL,OINFO,'Observation',ONAMES,OPTIONS)
%   creates a representation with specified observation names and options.
%   To create OPTIONS, use rlRepresentationOptions.
%
%   rep = rlRepresentation(ADMODEL,OINFO,AINFO,'Observation',ONAMES,'Action',ANames)
%   creates a representation for ADMODEL with default options and the
%   specified observation names. Additionally specify the network input
%   layer names, ANAMES, associated with each action specification as a
%   cell array. The names in ANAMES must be the same order as the action
%   specifications, AINFO.
%
%   rep = rlRepresentation(ADMODEL,OINFO,AINFO,'Observation',ONAMES,'Action',ANames,OPTIONS)
%   creates a representation with specified observation names, action
%   names, and options. To create OPTIONS, use rlRepresentationOptions.

%   Copyright 2018 The MathWorks, Inc.

    properties (Access = private)
        ExecutionSettings
        Assembler
        AnalyzedLayers
        Regularizer
        GradThresholder
        Solver
        Precision
        NetworkInfo
        NeedsStatefulTraining
        
        % network io names determined from the preBuildCheckNetwork
        NetworkInputNames
        NetworkOutputNames
        
        SaveData
    end
    
    methods
        %% Constructor
        function this = rlLayerRepresentation(rlFunction,varargin)
            this = this@rl.util.rlAbstractRepresentation(rlFunction,varargin{:});
            if isa(this.rlFunction,'nnet.cnn.layer.Layer')
                AssignedLayerNames = arrayfun(@(x) x.Name,this.rlFunction,'UniformOutput',false);
                ix = cellfun(@(x) isempty(x),AssignedLayerNames);                                
                AssignedLayerNames(ix) = matlab.lang.makeUniqueStrings(repmat({'LayerName'},sum(ix),1),AssignedLayerNames(~ix));                
                for ct = 1:length(ix)
                    this.rlFunction(ct).Name = AssignedLayerNames{ct};                    
                end
                
                this.rlFunction = layerGraph(this.rlFunction);
            end            
            
            this.LocalOptions = rl.util.rlLayerRepresentation.createOptimizerOptions(this.Options);
            
            % check the validity of the network before it is built
            this = preBuildCheckNetwork(this);
            
            if this.hasOutputLayers
                % the network is not complete if there is no output layer
                % wait for setLoss function to de
                this = buildNetwork(this);
                this = validateRepresentation(this);
            end
        end       
        
        %% Set loss
        function this = setLoss(this,CurrentLoss,varargin)
            CurrentLoss = validatestring(CurrentLoss,{'mse','cte','ppo-clipped'},'','Loss');
            switch CurrentLoss
                case {'cte','mse'}
                    DefaultEntropyLossWeight = 0;
                    DefaultClipFactor = 0;
                case 'ppo-clipped'
                    DefaultEntropyLossWeight = 0.01;
                    DefaultClipFactor = 0.2;
            end
            parser = inputParser;
            addParameter(parser,'EntropyLossWeight',DefaultEntropyLossWeight);
            addParameter(parser,'ClipFactor',DefaultClipFactor);
            parse(parser,varargin{:});
            result = parser.Results;
            EntropyLossWeight = result.EntropyLossWeight;
            ClipFactor = result.ClipFactor;
            validateattributes(EntropyLossWeight,{'numeric'},{'scalar','nonnegative','finite','real','<=',1},'','EntropyLossWeight');
            validateattributes(ClipFactor,{'numeric'},{'scalar','nonnegative','finite','real','<',1},'','ClipFactor');
            
            if this.hasOutputLayers
                % already loss layers are appended, 
                % don't set loss, just check it's compatible
                if this.Validated
                    layers  = this.adModel.OutputLayers;
                    switch CurrentLoss
                        case 'mse'
                            if ~any(cellfun(@(x) isa(x,'nnet.internal.cnn.layer.MeanSquaredError'),layers))
                                error(message('rl:agent:errRequiredMSELoss'));
                            end
                        case 'cte'
                            if ~any(cellfun(@(x) isa(x,'nnet.internal.cnn.layer.CrossEntropy'),layers))
                                ix = cellfun(@(x) isa(x,'nnet.internal.cnn.layer.CustomClassificationLayer'),layers);
                                % check whether attached loss is cte with entropy
                                if any(ix)
                                    % assume only classification layer with scalar EntropyWeightLoss
                                    ExternalLayer = layers{ix}.ExternalCustomLayer;
                                    if isa(ExternalLayer,'rl.layer.internal.ClassificationLayerWithEntropy')
                                        % set OutputLayers new entropy loss weight
                                        if EntropyLossWeight ~= ExternalLayer.EntropyLossWeight
                                            ExternalLayer.EntropyLossWeight = EntropyLossWeight;
                                            this.adModel.OutputLayers{ix} = nnet.internal.cnn.layer.CustomClassificationLayer( ...
                                                ExternalLayer, this.adModel.OutputLayers{ix}.LayerVerifier);
                                        end
                                    else
                                        error(message('rl:agent:errRequiredCTELoss'));
                                    end
                                else
                                    % error if attached loss is not cte or cte with entropy
                                    error(message('rl:agent:errRequiredCTELoss'));
                                end
                            else
                                if EntropyLossWeight
                                    % error if cte loss is already attached
                                    % and entropy loss is enabled
                                    error(message('rl:agent:errRemoveCTELoss'));
                                end
                            end
                        case 'ppo-clipped'
                            ix = cellfun(@(x) isa(x,'nnet.internal.cnn.layer.CustomClassificationLayer'),layers);
                            % check whether attached loss is cte with entropy
                            if any(ix)
                                % assume only classification layer with scalar EntropyWeightLoss
                                ExternalLayer = layers{ix}.ExternalCustomLayer;
                                if isa(ExternalLayer,'rl.layer.internal.ClipPPOLossLayer')
                                    % set OutputLayers new entropy loss weight
                                    if EntropyLossWeight ~= ExternalLayer.EntropyLossWeight
                                        ExternalLayer.EntropyLossWeight = EntropyLossWeight;
                                        this.adModel.OutputLayers{ix} = nnet.internal.cnn.layer.CustomClassificationLayer( ...
                                            ExternalLayer, this.adModel.OutputLayers{ix}.LayerVerifier);
                                    end
                                    % set OutputLayers new clip factor
                                    if ClipFactor ~= ExternalLayer.ClipFactor
                                        ExternalLayer.ClipFactor = ClipFactor;
                                        this.adModel.OutputLayers{ix} = nnet.internal.cnn.layer.CustomClassificationLayer( ...
                                            ExternalLayer, this.adModel.OutputLayers{ix}.LayerVerifier);
                                    end
                                else
                                    error(message('rl:agent:errRequiredPPOLoss'));
                                end
                            else
                                % error if attached loss is not ppo loss
                                error(message('rl:agent:errRequiredPPOLoss'));
                            end
                    end
                else
                    error(message('rl:agent:errLossDefined'));
                end
            else
                % get the digraph to the network
                dg = getDigraph(this);
                % then get the adjacency matrix
                adj = dg.adjacency';
                % find the output (sink) nodes by finding nodes that are
                % NOT inputs to any other node.
                snkidx = ~any(adj,1);
                % get the layer names
                Names = {this.rlFunction.Layers.Name};
                % get the output layer name
                OutputName = Names{snkidx};
                switch CurrentLoss
                    case 'mse'
                        % append mse loss
                        Names = matlab.lang.makeUniqueStrings([Names {'RepresentationLoss'}]);
                        CurrentLossName = Names{end};
                        CurrentInputNode = CurrentLossName;
                        CurrentLossLayer = regressionLayer('Name',CurrentLossName);
                    case {'cte', 'ppo-clipped'}
                        % append cte or ppo loss
                        if ~isa(this.ActionInfo,'rl.util.rlFiniteSetSpec')
                            error(message('rl:agent:errCTEDiscreteOutput'));
                        end
                        Names = matlab.lang.makeUniqueStrings([Names {'RepresentationSoftMax'} {'RepresentationLoss'}]);
                        CurrentLossName = Names{end};
                        CurrentSoftmaxName = Names{end-1};
                        CurrentInputNode = CurrentLossName;
                        CurrentLossLayer = [];
                        if ~isa(this.rlFunction.Layers(snkidx),'nnet.cnn.layer.SoftmaxLayer')
                            % if last layer is not softmax, append it
                            CurrentInputNode = CurrentSoftmaxName;
                            CurrentLossLayer = [CurrentLossLayer softmaxLayer('Name',CurrentSoftmaxName)];
                        end
                        switch CurrentLoss
                            case 'cte'
                                if EntropyLossWeight
                                    CurrentLossLayer = [CurrentLossLayer ...
                                        rl.layer.internal.ClassificationLayerWithEntropy('Name',CurrentLossName,'EntropyLossWeight',EntropyLossWeight)];
                                else
                                    CurrentLossLayer = [CurrentLossLayer classificationLayer('Name',CurrentLossName)];
                                end
                            case 'ppo-clipped'
                                CurrentLossLayer = [CurrentLossLayer ...
                                    rl.layer.internal.ClipPPOLossLayer('Name',CurrentLossName,'EntropyLossWeight',EntropyLossWeight,'ClipFactor',ClipFactor)];
                        end
                end
                tempRepresentation = this.rlFunction;
                tempRepresentation = addLayers(tempRepresentation,CurrentLossLayer);
                tempRepresentation = connectLayers(tempRepresentation,OutputName,CurrentInputNode);
                
                this.rlFunction = tempRepresentation;
                this = buildNetwork(this);
                this = validateRepresentation(this);
            end
        end
        
        %% Saving object
        function s = saveobj(this)
            % General properties
            s.Options = this.Options;
            
            if this.Validated
                % rep has loss function                                                
                if isActor(this)
                    % Remove all layers after the layer with action name
                    downstreamLayers = iAllDownstreamLayers(this.rlFunction, this.ActionNames{1});
                    s.rlFunction = removeLayers(this.rlFunction, downstreamLayers);
                    s.Validated = false;
                    % From version 2, s.Loss is a struct contains loss 
                    % name and loss parameters (if any). e.g 
                    % s.Loss.LossName = 'cte'; s.Loss.LossParam.EntropyLossWeight = 0.01;
                    % In version 1, s.Loss is loss name: e.g. s.Loss = 'mse'
                    s.Loss = getLoss(this);
                else
                    s.rlFunction = this.rlFunction;
                    s.Validated = true;
                    s.Loss = getLoss(this);
                end
                % Serialize the network
                s.adModelExternal = this.getExternalNetwork();
            else
                % rep no loss, just save rlFunction
                s.rlFunction = this.rlFunction;
                s.Validated = false;
                s.Loss = []; % we don't know the loss
                s.adModelExternal = [];
            end
                        
            s.ObservationInfo = this.ObservationInfo;
            s.ActionInfo = this.ActionInfo;
            s.ObservationNames = this.ObservationNames;
            s.ActionNames = this.ActionNames;
            s.SaveData = this.SaveData;
                                            
            % save representation with UseDevice set to 'cpu'
            s.SaveData.UseDevice = this.Options.UseDevice;
            s.Options.UseDevice = 'cpu';
            
            % version indicator
            s.Version = 2;
        end

        %% Get external network
        function mdl = getModel(this)
            % Return the computation model as a SeriesNetwork or DAGNetwork
            % neural network object.
            
            if this.hasOutputLayers
                mdl = getExternalNetwork(this);
            else
                error(message('rl:agent:errLayerGetModelLossRequired'))
            end
        end
    end
    methods (Hidden)      
        function Parameters = getParameters(this)
            Parameters = this.adModel.LearnableParameters;
        end
        
        %% Compute gradients of representation output wrt learnable parameters or inputs
        function [gradVal,paramVal,stateVal] = gradient(this,dy,dx,inputValues,varargin)
            % Required: (dOutput/dParameter,dOutput/dAction,dOutput/dObservation,dLoss/dParameter)
            if nargin<6
                initialGrad = false;
            else
                initialGrad = true;
            end
            
            % get output and input nodes            
            dy = lower(dy);
            dx = lower(dx);
            
            Data = mapInputValues(this,inputValues);
            Data = cellfun(@(x) processData(this,x),Data,'UniformOutput',false);
            
            % compute gradients
            if strcmpi(dx,'observation') || strcmpi(dx,'action')
                nBatch = size(Data{1},4);
                [~, ~, stateVal, gradVals] = this.adModel.computeGradientsForTraining( ...
                         Data, [], false,ones([this.adModel.OutputSizes{1} nBatch],'single'));
                if strcmpi(dx,'observation')
                    gradVal = gradVals(this.ObservationIndex);
                elseif strcmpi(dx,'action') && strcmpi(this.ActionIndex{1},'input')
                    gradVal = gradVals(this.ActionIndex{2});
                else
                    gradVal = [];
                end                     
                paramVal = [];
            elseif strcmpi(dx,'parameter')
                if initialGrad
                    % reshape initial gradients to match actor's output size
                    initialGradients = cellfun(@(x,y) reshape(x,[y size(Data{1},4)]),...
                        varargin(2),this.adModel.OutputSizes,'UniformOutput',false);
                    initialGradients = cellfun(@(x) processData(this,x),initialGradients,'UniformOutput',false);
                    [gradVal,~,stateVal] = this.adModel.computeGradientsForTraining(Data,[],false,initialGradients);
                    paramVal = [];
                else
                    LossVariable = varargin{1};
                    switch numel(LossVariable)
                        case 1
                            LossVariable = processData(this,LossVariable{1});
                            LossVariable = convertTargetData(this,LossVariable);
                        otherwise
                            LossVariable = {LossVariable};
                    end
                    [gradVal,~,stateVal] = this.adModel.computeGradientsForTraining(Data,LossVariable,false);
                    paramVal = [];
                end
            end            
        end
        
        %% Step learnable parameters with input and target data
        function this = fit(this,InputData,LossVariable)
            
            % only compute gradients and update params if LR > 0
            if this.Options.LearnRate > 0
                [gradVal, ~, stateVal] = gradient(this,'loss','parameter',InputData,LossVariable);
                this = step(this,gradVal,stateVal);
            end
        end
        
        %% Step learnable parameters with gradient information or loss
        function this = step(this,gradient,state)
            % only update params if LR > 0
            if this.Options.LearnRate > 0
                gradient = regularizeGradients(this,gradient);
                gradient = thresholdGradients(this,gradient);
                this = applyGradients(this,gradient,state);
            end
        end
        
        %% Regularize gradient
        function gradient = regularizeGradients(this,gradient)
            gradient = this.Regularizer.regularizeGradients(gradient,this.adModel.LearnableParameters);
        end
        
        %% Threshold gradient
        function gradient = thresholdGradients(this,gradient)
            gradient = thresholdGradients(this.GradThresholder,gradient);            
        end
        
        %% Apply gradient
        function this = applyGradients(this,gradient,state)
            velocity = this.Solver.calculateUpdate(gradient,this.LocalOptions.InitialLearnRate);
            this.adModel = this.adModel.updateLearnableParameters(velocity);
            this.adModel = this.adModel.updateNetworkState(state);            
        end
        
        %% Copy representation
        function Representation2 = copy(Representation1)
            % Representation is value object
            Representation2 = Representation1;
        end
        
        %% Return rlFunction (it may be different than input argument)
        function net = getExternalNetwork(this)
           TempadModel = this.adModel.prepareNetworkForPrediction();
           TempadModel = TempadModel.setupNetworkForHostPrediction();
           net = this.Assembler.createExternalNetwork(TempadModel,this.AnalyzedLayers);  
        end
        function rlFunc = qeGetRLFunction(this)
           rlFunc = this.rlFunction;
        end
        function rlFunc = qeGetADModel(this)
           rlFunc = this.adModel;
        end
        function outputSize = getLayerOutputSizeFromName(this, name)
            % THIS CODE IS BORROWED FROM DAGNetwork since the change is
            % submitted after CF and we need to test in our code in
            % Bcontrol
            layerNames = cellfun(@(x)x.Name, this.adModel.SortedLayers, 'UniformOutput', false);
            sortedLayerIndex = find(strcmp(name, layerNames));
            outputSize = this.adModel.inferOutputSizesGivenInputSizes(this.adModel.InputSizes, sortedLayerIndex); %#ok<FNDSB>
        end
        
        %% Get loss information
        function Loss = getLoss(this)
            % Traverse the network to collect loss information
            % Return a struct which contains loss name and loss parameters
            % e.g Loss.LossName = 'cte'; Loss.EntropyLossWeight = 0.01;
            
            if any(cellfun(@(x) isa(x,'nnet.internal.cnn.layer.MeanSquaredError'),this.adModel.OutputLayers))
                LossName = 'mse';
                LossParam = [];
            elseif any(cellfun(@(x) isa(x,'nnet.internal.cnn.layer.CrossEntropy'),this.adModel.OutputLayers))
                LossName = 'cte';
                LossParam.EntropyLossWeight = 0;
            end
            ix = cellfun(@(x) isa(x,'nnet.internal.cnn.layer.CustomClassificationLayer'),this.adModel.OutputLayers);
            if any(ix)
                % assume only one custom loss layer per representation
                ExternalLayer = this.adModel.OutputLayers{ix}.ExternalCustomLayer;
                switch class(ExternalLayer)
                    case 'rl.layer.internal.ClassificationLayerWithEntropy'
                        LossName = 'cte';
                        LossParam.EntropyLossWeight = ExternalLayer.EntropyLossWeight;
                    case 'rl.layer.internal.ClipPPOLossLayer'
                        LossName = 'ppo-clipped';
                        LossParam.EntropyLossWeight = ExternalLayer.EntropyLossWeight;
                        LossParam.ClipFactor = ExternalLayer.ClipFactor;
                end
            end
            Loss.LossName = LossName;
            Loss.LossParam = LossParam;
        end
    end
    methods (Access = protected)
        function Value = evaluate_(this,Data)
            % Map to network inputs
            Data = mapInputValues(this,Data);
            Data = cellfun(@(x) processData(this,x),Data,'UniformOutput',false);
            Value = predict(this.adModel,Data);
            Value = convertOutputData(this,Value);
            % TODO: there should be mapping for output values if action is at the
            % output
        end
        function Value = getLearnableParameterValues_(this)
            params = getParameters(this);
            for ct=1:numel(params)
                Value{1,ct} = params(ct).Value;
            end            
        end
        function this = setLearnableParameterValues_(this,Value)
            this.adModel = setLearnableParameterValues(this.adModel,Value);            
        end
        function sz = getSize_(this,identifier)
            switch identifier
                case 'observation'
                    Inputs = this.adModel.InputSizes;
                    ix = this.ObservationIndex;
                    for ct=1:length(ix)
                        sz{1,ct} = Inputs{ix(ct)};
                    end
                case 'action'
                    sz = {};
                    if ~isempty(this.ActionNames)
                        ix = this.ActionIndex{2};
                        if strcmpi(this.ActionIndex{1},'input')
                            sz = this.adModel.InputSizes(ix);
                        else
                            sz = this.getOutputActionSizes;
                        end                     
                    end
                case 'output'
                    Outputs = this.adModel.OutputSizes;
                    for ct=1:numel(Outputs)
                        sz{1,ct} = [Outputs{ct}(end) 1];
                    end
                case 'input'
                    sz = this.adModel.InputSizes;
                    % add input order: IMP
                case 'parameter'
                    sz = cellfun(@(x) size(x),getLearnableParameterValues(this),'UniformOutput',false);
            end
        end
        function validateLearnableParameterSize(this,Value)
            if ~isa(Value,'cell')
                error(message('rl:agent:errLearnableParameters'));
            end
            params = getParameters(this);
            paramSize = arrayfun(@(x) size(x.Value),params,'UniformOutput',false);
            valSize = cellfun(@(x) size(x),Value,'UniformOutput',false);
            if length(paramSize)~=length(valSize)
                error(message('rl:agent:errParameterSize'));                
            end
            if ~all(cellfun(@(x,y) isequal(x,y),paramSize,valSize))
                error(message('rl:agent:errParameterSize'));
            end            
        end
        function this = validateRepresentation(this)

            InputNames = reshape(this.NetworkInputNames,1,[]);    
            OutputNames = reshape(this.NetworkOutputNames,1,[]);
            
            %% find observation name, index and get its size
            ObservationNames = this.ObservationNames;
            [tmp,~,ix] = intersect(ObservationNames,InputNames,'stable');
            ix = reshape(ix,1,[]);
            if isempty(tmp)
                error(message('rl:agent:errObservationNames'));
            end
            this.ObservationIndex = ix;
            ObservationSizes = this.adModel.InputSizes(ix);
            this.InputIndex = ix;
            
            %% find action name, index and get its size
            ActionNames = this.ActionNames;
            ActionSizes = [];
            if ~isempty(ActionNames)
                if isempty(intersect(ActionNames,union(InputNames,OutputNames),'stable'))
                    error(message('rl:agent:errActionNames'));
                end
                
                [tmp,~,ix] = intersect(ActionNames,InputNames,'stable');
                ix = reshape(ix,1,[]);
                if ~isempty(tmp)
                    this.ActionIndex = {'input',ix};
                    ActionSizes = reshape(this.adModel.InputSizes(ix),1,[]);
                    this.InputIndex = [this.InputIndex this.ActionIndex{2}];
                end
                
                [tmp,~,ix] = intersect(ActionNames,OutputNames,'stable');
                ix = reshape(ix,1,[]);
                if ~isempty(tmp)
                    this.ActionIndex = {'output',ix};
                    ActionSizes = getOutputActionSizes(this);
                end                                
            end
            
            %% validate number of imageInputLayers is equal to observation/action
            numInputs = this.adModel.NumInputLayers;                        
            numTotalInfo = numel(this.ObservationInfo);
            if ~isActor(this)
                numTotalInfo = numTotalInfo + numel(this.ActionInfo);
            end
            if numInputs ~= numTotalInfo
                error(message('rl:agent:errIncompatibleNumberOfImageInputLayer',numInputs,numTotalInfo));
            end                        
            
            validateRepresentationSizes(this,ObservationSizes,ActionSizes);            
            this.Validated = true;
            
            %% usample and validate network by forward pass
            try
                % not to modify the random seed, get state and set it back
                % after usample.
                inputVals = usampleInputs(this);
                val = evaluate(this,inputVals);
                % check the actor output case (not critic output since it
                % may be restrictive and checked in agent level)
                if isActor(this)
                    outSize = getSize(this,'output');
                    if ~iscell(val)
                        val = {val};
                    end
                    valSize = cellfun(@(x) size(x),val,'UniformOutput',false);
                    if ~all(cellfun(@(x,y) rl.util.isArrayEqual(x,y),valSize,outSize))
                        this.Validated = false;
                        error(message('rl:agent:errIncompatibleAction2'))
                    end
                end                
            catch ex
                this.Validated = false;
                me = MException(message('rl:agent:errRepresentationEvaluation'));
                throw(addCause(me,ex))
            end
        end         
        function this = preBuildCheckNetwork(this)
            % do some checks on the network before it is "built"
            
            % get the digraph
            [dg,layerNodes] = getDigraph(this);
            
            %% make sure the graph is DAG
            if ~isdag(dg)
                % TODO MSG
                error(message('rl:agent:errNotDAG'));
            end
            
            
            %% make sure all layers are connected (single bin)
            bins = conncomp(dg,'Type','weak');
            badbins = bins > 1;
            if any(badbins)
                badlayers = arrayfun(@(x) x.Name,layerNodes(badbins),'UniformOutput',false);
                balayersstr = sprintf('"%s", ',badlayers{:});
                balayersstr(end-1:end) = '';
                balayersstr = ['[',balayersstr,']'];
                error(message('rl:agent:errDisconnectedLayers',balayersstr));
            end
            
            %% check for any unsupported layers
            if any(arrayfun(@(x) ...
                    isa(x,'nnet.cnn.layer.SequenceInputLayer') || ...
                    isa(x,'nnet.cnn.layer.LSTMLayer') || ...
                    isa(x,'nnet.cnn.layer.BatchNormalizationLayer') || ...
                    isa(x,'nnet.cnn.layer.BiLSTMLayer'),layerNodes))
                error(message('rl:agent:errInvalidLayers'));
            end
            
            %% make sure that Normalization is set to None in imageInputLayer
            inputLayer = layerNodes(arrayfun(@(x) isa(x,'nnet.cnn.layer.ImageInputLayer'),layerNodes));
            if any(arrayfun(@(x) ~strcmpi(x.Normalization,'none'),inputLayer))
                error(message('rl:agent:errImageInputLayerNormalizationNone'));
            end
            
            %% make sure the IOs are defined
            % get the adjacency matrix
            adj = dg.adjacency';
            % find the output (sink) nodes by finding nodes that are
            % NOT inputs to any other node.
            yidx = find(~any(adj,1));
            for i = 1:numel(yidx)
                yi = yidx(i);
                outLayer = layerNodes(yi);
                if isa(outLayer,'nnet.internal.cnn.layer.ClassificationLayer') || isa(outLayer,'nnet.cnn.layer.ClassificationOutputLayer') || isa(outLayer, 'rl.layer.internal.ClassificationLayerWithEntropy') || isa(outLayer, 'rl.layer.internal.ClipPPOLossLayer')
                    % CTE: One layer back due to softmax
                    ix = predecessors(dg,yi);
                    layer = layerNodes(ix);
                    if ~ ( isa(layer,'nnet.internal.cnn.layer.SoftmaxLayer') || isa(layer,'nnet.cnn.layer.SoftmaxLayer') )
                        error(message('rl:agent:errMissingSoftmaxLayer'));
                    end                    
                    yidx(i) = ix;
                elseif isa(outLayer,'nnet.internal.cnn.layer.SoftmaxLayer') || isa(outLayer,'nnet.cnn.layer.SoftmaxLayer')
                    yidx(i) = yi;
                elseif isa(outLayer,'nnet.internal.cnn.layer.RegressionLayer') || isa(outLayer,'nnet.cnn.layer.RegressionOutputLayer')
                    % MSE
                    yidx(i) = predecessors(dg,yi);
                end
            end
            % find input (src) nodes by finding nodes that have NO incoming
            % connections
            uidx = ~any(adj,2)';
            
            % get the io names
            fcn = @(x) x.Name;
            yNames = arrayfun(fcn,layerNodes(yidx),'UniformOutput',false);
            uNames = arrayfun(fcn,layerNodes(uidx),'UniformOutput',false);
            
            % find observation name, index and get its size
            ObservationNames = this.ObservationNames;
            if isempty(intersect(ObservationNames,uNames,'stable'))
                error(message('rl:agent:errObservationNames'));
            end
            
            % find action name, index and get its size
            ActionNames = this.ActionNames;
            if ~isempty(ActionNames)
                if isempty(intersect(ActionNames,union(uNames,yNames),'stable'))
                    error(message('rl:agent:errActionNames'));
                end                              
            end 
            
            this.NetworkOutputNames = yNames;
            this.NetworkInputNames  = uNames;
        end
        function Data = convertTargetData(this,Data)
            if iscell(Data)
                for ct=1:numel(Data)
                    sz = size(Data{ct});
                    Data{ct} = reshape(Data{ct},[1 1 sz(1) sz(end)]);
                end
            else
                sz = size(Data);
                Data = reshape(Data,[1 1 sz(1) sz(end)]);
            end
        end
        function Data = convertOutputData(this,Data)
            if iscell(Data)
                for ct=1:numel(Data)
                    sz = size(Data{ct});
                    sz = [sz ones(1,4-length(sz))];
                    Data{ct} = reshape(Data{ct},[sz(3) 1 1 sz(end)]);
                end
            else
                sz = size(Data);
                sz = [sz ones(1,4-length(sz))];
                Data = reshape(Data,[sz(3) 1 1 sz(end)]);
            end
            
            if iscell(Data)
                Data = Data{1};                
            end
            
            if isGPU(this)
                Data = gather(Data);
            end
        end
        function [dg,layerNodes] = getDigraph(this)
            % get the digraph to the network
            dg = extractPrivateDirectedGraph(this.rlFunction);
            layerNodes = dg.Nodes{:,1};
        end
        function this = updateOptions(this,opt)
            % convert representation option to local options
            this.LocalOptions = rl.util.rlLayerRepresentation.createOptimizerOptions(opt);                        
            if this.Validated                
                copt = this.Options;                               
                this.ExecutionSettings = nnet.internal.cnn.assembler.setupExecutionEnvironment(this.LocalOptions, this.NetworkInfo.IsRNN, [], this.Precision );
                this.Regularizer = nnet.internal.cnn.regularizer.RegularizerFactory.create( ...
                                   'l2', this.adModel.LearnableParameters, this.Precision, this.LocalOptions );
                this.Solver = nnet.internal.cnn.solver.SolverFactory.create( ...
                              this.adModel.LearnableParameters, this.Precision, this.LocalOptions );
                
                % check whether gradient threshold options are changed
                if ~( isequal(copt.GradientThreshold,opt.GradientThreshold) && ...
                        isequal(copt.GradientThresholdMethod,opt.GradientThresholdMethod) )
                    gradThresOptions = struct('Method',opt.GradientThresholdMethod, ...
                        'Threshold',opt.GradientThreshold);
                    this.GradThresholder = nnet.internal.cnn.GradientThresholder(gradThresOptions);
                end 
                
                this = prepareADModelForTraining(this);
            end
        end
        function argStruct = generateEvaluateFunction_(this,argStruct)
            % generate the evaluate function for the layer representation
            outputstr = argStruct.EvaluateOutputString;
            inputstr  = argStruct.EvaluateInputString;
            policyname = argStruct.PolicyName;
            matfilename = argStruct.MATFileName;
            localfcnstr = argStruct.LocalFunctionString;
            evalfcnname = argStruct.EvaluateFunctionName;

            %% Extract vars and make sure the DAG can be constructed
            oinfo = this.ObservationInfo;
            ainfo = this.ActionInfo;
            
            if ~isActor(this) || ~isscalar(ainfo) || ~isscalar(oinfo)
                error(message('rl:general:RepLayerCodeGenMIMONotSupported'));
            end
            
            %% Extract the DAG
            dag = getExternalNetwork(this); %#ok<NASGU>
            
            %% Save the network to disk
            eval(sprintf('%s = dag;',policyname));
            save(matfilename,policyname);
            
            %% generate the local function
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
            
            %% attach it to the arg struct
            if isempty(localfcnstr)
                localfcnstr = localevalfcnstr;
            else
                localfcnstr = sprintf('%s\n%s',localfcnstr,localevalfcnstr);
            end
            argStruct.LocalFunctionString = localfcnstr;
        end
    end
    methods (Access = private)
        function [internalNetwork,executionSettings,assembler,analyzedLayers,networkInfo] = createInternalNeuralNetwork(this)
            lgraph = this.rlFunction;
            
            %% Check network for SISO & MIMO 
            isaDAG = isa(lgraph,'nnet.cnn.LayerGraph');
            analyzedLayers = nnet.internal.cnn.analyzer.NetworkAnalyzer(lgraph);
            % nnet.internal.halfbaked.mimo.nnet.internal.cnn.analyzer.constraints.Architecture()
            analyzedLayers.applyConstraints([
                rl.util.Architecture()
                nnet.internal.cnn.analyzer.constraints.Connections()
                nnet.internal.cnn.analyzer.constraints.CustomLayers()
                nnet.internal.cnn.analyzer.constraints.LSTM()
                nnet.internal.cnn.analyzer.constraints.Names()
                nnet.internal.cnn.analyzer.constraints.Propagation()]);
            analyzedLayers.throwIssuesIfAny()
            internalLayers = analyzedLayers.InternalLayers;
            networkInfo = nnet.internal.cnn.util.ComputeNetworkInfo(isaDAG,internalLayers);
            
            % Set up and validate training
            executionSettings = nnet.internal.cnn.assembler.setupExecutionEnvironment(...
                this.LocalOptions, networkInfo.IsRNN, [], this.Precision );
            
            % Assemble internal network
            strategy = nnet.internal.cnn.assembler.NetworkAssemblerStrategyFactory.createStrategy(~networkInfo.IsDAG); % try with removing it
            assembler = nnet.internal.cnn.assembler.TrainingNetworkAssembler(strategy);            
            
            internalNetwork = assembler.assemble(analyzedLayers, executionSettings);
        end
        function this = buildNetwork(this)
            % Set desired precision (not move to options)
            this.Precision = nnet.internal.cnn.util.Precision('single');            
            
            [this.adModel, this.ExecutionSettings, ...
             this.Assembler, this.AnalyzedLayers,this.NetworkInfo] = createInternalNeuralNetwork(this);
                                   
            % TODO: Check MIMO loss function case
            % Assume that this network has a single classification output layer. Store
            % required classification meta-data in that layer.
            for i = 1:this.adModel.NumOutputs
                if isa(this.adModel.OutputLayers{i}, 'nnet.internal.cnn.layer.CrossEntropy') || ...
                   isa(this.adModel.OutputLayers{i}, 'nnet.internal.cnn.layer.CustomClassificationLayer')
                    if isempty(this.ActionInfo) || ~isa(this.ActionInfo,'rl.util.rlFiniteSetSpec')
                        error(message('rl:agent:errCTEDiscreteOutput'));
                    end
                    
                    ActionValues = reshape(this.ActionInfo.Elements,[],1);
                    % Single channel, multi actions case:
                    % ActionValues is a cell array, cast to string before
                    % convert to meta-data.
                    if iscell(ActionValues)
                        ActionValues = cellfun(@mat2str,ActionValues,'UniformOutput',false);
                    end
                    this.adModel.OutputLayers{i}.Categories = ...
                        categorical(ActionValues, ActionValues);
                end
            end
            
            % Create L2 regularizer
            this.Regularizer = nnet.internal.cnn.regularizer.RegularizerFactory.create( ...
                'l2', this.adModel.LearnableParameters, this.Precision, this.LocalOptions );
            
            % Create solver
            this.Solver = nnet.internal.cnn.solver.SolverFactory.create( ...
                this.adModel.LearnableParameters, this.Precision, this.LocalOptions );
            
            % Create gradient thresholder
            gradientThresholdOptions = struct('Method', this.LocalOptions.GradientThresholdMethod,...
                'Threshold', this.LocalOptions.GradientThreshold);
            this.GradThresholder = nnet.internal.cnn.GradientThresholder(gradientThresholdOptions);
            
            % Book keeping
            this.NeedsStatefulTraining = nnet.internal.cnn.util.isStatefulLayer( this.adModel.Layers );            
        end
        function flag = hasOutputLayers(this)
            ix = arrayfun(@(x) isa(x,'nnet.cnn.layer.RegressionOutputLayer') || ...
                               isa(x,'nnet.cnn.layer.ClassificationOutputLayer') || ...
                               isa(x,'rl.layer.internal.ClassificationLayerWithEntropy') || ...
                               isa(x,'rl.layer.internal.ClipPPOLossLayer'),this.rlFunction.Layers);
            flag = any(ix);
        end
        function sz = getOutputActionSizes(this)
            ActionNames = this.ActionNames;
            for ct=1:numel(ActionNames)
                 tmp = getLayerOutputSizeFromName(this,ActionNames{ct});
                 sz{1,ct} = [tmp{1}(end) 1];
            end
        end
        function Data = processData(this,Data)
            if isGPU(this)
                Data = gpuArray(this.Precision.cast(Data));
            else
                Data = this.Precision.cast(gather(Data));
            end
        end
        function flag = isGPU(this)
            flag = ismember(this.ExecutionSettings.executionEnvironment, {'gpu'} );
        end
        function this = prepareADModelForTraining(this)
            % Convert learnable params to training mode with correct
            % excecution environment (e.g. convert to gpuArray for GPU usage)
            this.ExecutionSettings.executionEnvironment = this.LocalOptions.ExecutionEnvironment;
            this.ExecutionSettings.useParallel = false;
            this.adModel = this.adModel.prepareNetworkForTraining(this.ExecutionSettings);
        end
    end
    methods (Static,Hidden)
        function newopt = createOptimizerOptions(opt)
            % supported optimizers
            optimizer = validatestring(opt.Optimizer,{'adam','sgdm','rmsprop'},'','Optimizer');
            % if the learn rate is set to zero, fake it for DLTs sake. fit
            % and step will bypass learning in this case
            if opt.LearnRate == 0
                opt.LearnRate = eps;
            end
            switch optimizer
                case 'adam'
                    newopt = trainingOptions(opt.Optimizer,...
                        'InitialLearnRate',opt.LearnRate, ...
                        'GradientDecayFactor',opt.OptimizerParameters.GradientDecayFactor, ...
                        'SquaredGradientDecayFactor',opt.OptimizerParameters.SquaredGradientDecayFactor, ...
                        'Epsilon',opt.OptimizerParameters.Epsilon, ...
                        'GradientThresholdMethod',opt.GradientThresholdMethod, ...
                        'GradientThreshold',opt.GradientThreshold, ...
                        'L2Regularization',opt.L2RegularizationFactor, ...
                        'ExecutionEnvironment',char(opt.UseDevice));
                case 'sgdm'
                    newopt = trainingOptions(opt.Optimizer,...
                        'InitialLearnRate',opt.LearnRate, ...
                        'Momentum',opt.OptimizerParameters.Momentum, ...
                        'GradientThresholdMethod',opt.GradientThresholdMethod, ...
                        'GradientThreshold',opt.GradientThreshold, ...
                        'L2Regularization',opt.L2RegularizationFactor, ...
                        'ExecutionEnvironment',char(opt.UseDevice));
                case 'rmsprop'
                    newopt = trainingOptions(opt.Optimizer,...
                        'InitialLearnRate',opt.LearnRate, ...                        
                        'SquaredGradientDecayFactor',opt.OptimizerParameters.SquaredGradientDecayFactor, ...
                        'Epsilon',opt.OptimizerParameters.Epsilon, ...
                        'GradientThresholdMethod',opt.GradientThresholdMethod, ...
                        'GradientThreshold',opt.GradientThreshold, ...
                        'L2Regularization',opt.L2RegularizationFactor, ...
                        'ExecutionEnvironment',char(opt.UseDevice));
            end
        end 
        function validateRepresentationType(rep)
            if ~isa(rep,'rl.util.rlLayerRepresentation')
                error(message('rl:agent:errInvalidRepresentation'));
            end
        end
    end
    methods (Static)
        function this = loadobj(s)
            if ~isa(s,'rl.util.rlLayerRepresentation')
                s = iUpdatePreviousVersionSaveStruct(s);
                % General properties
                inputArgs = { ...
                    s.rlFunction, ...
                    s.ObservationInfo};
                if ~isempty(s.ActionInfo)
                    inputArgs{end+1} = s.ActionInfo;
                end
                if ~isempty(s.ObservationNames)
                    inputArgs{end+1} = 'Observation';
                    if iscell(s.ObservationNames) && (numel(s.ObservationNames) == 1)
                        inputArgs{end+1} = s.ObservationNames{1};
                    else
                        inputArgs{end+1} = s.ObservationNames;
                    end
                end
                if ~isempty(s.ActionNames)
                    inputArgs{end+1} = 'Action';
                    if iscell(s.ActionNames) && (numel(s.ActionNames) == 1)
                        inputArgs{end+1} = s.ActionNames{1};
                    else
                        inputArgs{end+1} = s.ActionNames;
                    end
                end
                inputArgs{end+1} = s.Options;
                this = rl.util.rlLayerRepresentation(inputArgs{:});                
                
                if this.Validated
                    % Load the network
                    this.adModel = iCreateInternalNetwork(s.adModelExternal);
                else
                    if ~isempty(s.Loss)
                        LossName  = s.Loss.LossName;
                        LossParam = s.Loss.LossParam;
                        % Compile actor with removed layers
                        switch s.Loss.LossName
                            case 'mse'
                                this = setLoss(this,LossName);
                            case 'cte'
                                this = setLoss(this,LossName,'EntropyLossWeight',LossParam.EntropyLossWeight);
                            case 'ppo-clipped'
                                this = setLoss(this,LossName,...
                                    'EntropyLossWeight',LossParam.EntropyLossWeight,...
                                    'ClipFactor',LossParam.ClipFactor);
                        end
                        
                        % Update adModel to ensure consistent with rlFunction
                        this.adModel = iCreateInternalNetwork(s.adModelExternal);
                    end
                end
                
                if ~isempty(s.Loss)
                    % Convert learnable params to training mode
                    % Not applicable for representation withwout loss
                    this = prepareADModelForTraining(this);
                end
                
                % revert to gpuArray if UseDevice is 'gpu'
                if isfield(this,'SaveData') || isprop(this,'SaveData')
                    if ~isempty(this.SaveData) || ~isempty(s.SaveData)
                        this.Options.UseDevice = s.SaveData.UseDevice;
                    end
                end
            else
                this = s;
                % Convert learnable params to training mode
                this = prepareADModelForTraining(this);
                % revert to gpuArray if UseDevice is 'gpu'
                if isfield(this,'SaveData') || isprop(this,'SaveData')
                    if ~isempty(this.SaveData)
                        this.Options.UseDevice = s.SaveData.UseDevice;
                    end
                end
            end
        end
    end    
end

function internalNetwork = iCreateInternalNetwork(externalNetwork)
if isa(externalNetwork, 'DAGNetwork')
    lgraph = layerGraph(externalNetwork);
else
    lgraph = nnet.cnn.LayerGraph(externalNetwork.Layers);
end
internalLayerGraph = iExternalToInternalLayerGraph( lgraph );
[internalLayerGraph, topologicalOrder] = internalLayerGraph.toposort(); 
internalNetwork = nnet.internal.cnn.DAGNetwork(internalLayerGraph, topologicalOrder);
end

function internalLayerGraph = iExternalToInternalLayerGraph( externalLayerGraph )
internalLayers = iGetInternalLayers( externalLayerGraph.Layers );
hiddenConnections = externalLayerGraph.HiddenConnections;
internalConnections = iHiddenToInternalConnections( hiddenConnections );
internalLayerGraph = nnet.internal.cnn.LayerGraph(internalLayers, internalConnections);
end

function internalConnections = iHiddenToInternalConnections( hiddenConnections )
internalConnections = nnet.internal.cnn.util.hiddenToInternalConnections( hiddenConnections );
end

function internalLayers = iGetInternalLayers( layers )
internalLayers = nnet.internal.cnn.layer.util.ExternalInternalConverter.getInternalLayers( layers );
end

function downstreamLayers = iAllDownstreamLayers(lg, layerName)

% Find all connections that have this layer as a source
layerNameWithSlash = [layerName '/'];
layerAsSourceIndices = strcmp(lg.Connections.Source, layerName) ...
    | strncmp(lg.Connections.Source, layerNameWithSlash, length(layerNameWithSlash));

% Find corresponding destinations
rawDestinations = lg.Connections.Destination(layerAsSourceIndices);

% Remove port information
destinations = extractBefore(rawDestinations, '/');
emptyDestinations = cellfun(@(x)isempty(x), destinations);
destinations(emptyDestinations) = rawDestinations(emptyDestinations);

% Remove duplicates
destinations = unique(destinations);

downstreamLayers = destinations;

if ~isempty(destinations)
    for i = 1:numel(destinations)
        downstreamLayers = [downstreamLayers; ...
            iAllDownstreamLayers(lg, destinations{i})]; %#ok<AGROW>
    end
end

% Remove duplicates
downstreamLayers = unique(downstreamLayers);

end

function LossParam = iGetLossParamFromExternalNet(ExternalNet)
% Only execute when load representation from version 1.
% Get loss parameters (use in setLoss) from external network
% Regression output will likely not hit, include for completion

if any(arrayfun(@(x) isa(x,'nnet.cnn.layer.RegressionOutputLayer'),ExternalNet.Layers))
    LossParam = [];
elseif any(arrayfun(@(x) isa(x,'nnet.cnn.layer.ClassificationOutputLayer'),ExternalNet.Layers))
    LossParam.EntropyLossWeight = 0;
else
    ix = arrayfun(@(x) isa(x,'rl.layer.internal.ClassificationLayerWithEntropy'),ExternalNet.Layers);
    LossLayer = ExternalNet.Layers(ix);
    LossParam.EntropyLossWeight = LossLayer.EntropyLossWeight;
end

end

function s = iUpdatePreviousVersionSaveStruct(s)
% Update previous version saveobj struct, before loading operation 

if ~isfield(s,'Version') && ~isempty(s.Loss)
    % Version 1 does not have Version field
    % From version 2, s.Loss is a struct contains loss 
    % name and loss parameters (if any). e.g.
    % s.Loss.LossName = 'cte'; 
    % s.Loss.LossParams.EntropyLossWeight = 0.01;
    % Recover loss params from external DLT network
    LossName = s.Loss;
    LossParam = iGetLossParamFromExternalNet(s.adModelExternal);
    s.Loss = struct('LossName',LossName,'LossParam',LossParam);
    s.Version = 2;
end
end