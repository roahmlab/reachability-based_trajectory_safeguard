classdef rlDLNetworkRepresentation < rl.util.rlAbstractRepresentation
    %rlDLNetworkRepresentation RL dlnetwork representation
    
    %   Copyright 2019 The MathWorks, Inc.
    
    properties (Access = private, Transient)
        % Solver object (to be replaced with DLT external solver)
        % This object is built in set.Options in rlAbstractRepresentation
        % Thus, it will be rebuild automatically when load
        RLSolver
        % Function handle of the network loss
        LossFcn
    end
    
    properties (Access = private)
        % Loss details (struct with loss name, parameters)
        Loss
        
        % Cache network IO size (cell arrays) from preBuildCheckNetwork()
        % Network input size
        InputSize
        % Network output size
        OutputSize
        % Representation observation size
        ObservationSizes
        % Representation action size
        ActionSizes
        % Names array of input layers
        NetworkInputNames
        % Names array of output layers
        NetworkOutputNames
        % Network output layers
        NetworkOutputLayers
    end
    
    methods
        %% Constructor
        function this = rlDLNetworkRepresentation(rlFunction,varargin)
            % rlDLNetworkRepresentation creates RL representation for RL 
            % agents using a neural network defined as a dlnetwork
            % from Deep Learning Toolbox

            % solver options and solver object are created in set.Options 
            % of rlAbstractRepresentation constructor
            this = this@rl.util.rlAbstractRepresentation(rlFunction,varargin{:});
            
            % check the validity of the network before it is built
            this = preBuildCheckNetwork(this);
            
            % check the validity of the network by running a dummy predict
            % pass from observation and action info
            this = validateRepresentation(this);
        end
        
        %% Set loss
        function this = setLoss(this,CurrentLoss,varargin)
            % Set the loss function for the network
            
            CurrentLoss = validatestring(CurrentLoss,{'mse','cte','ppo-clipped'},'','Loss');
            % Input parser
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
            
            % Register network loss
            % LossVariable: cell array contains any necessary variables to
            % compute the loss.
            % e.g. For MSE loss, LossVariable contains target of the prediction
            % REVISIT: consider package loss constants (e.g. EntropyLossWeight)
            % to ConstantVariable (cell array) for consistency with LossVariable
            
            % Add softmax to dlnetwork if use classification loss
            if ismember(CurrentLoss,{'cte','ppo-clipped'})
                if ~isa(this.ActionInfo,'rl.util.rlFiniteSetSpec')
                    error(message('rl:agent:errCTEDiscreteOutput'));
                end
                if isempty(this.LossFcn)
                    % Only attach softmaxLayer if Loss is not attached
                    % REVISIT: currenly assume single output/ single loss
                    if ~any(arrayfun(@(x) isa(x,'nnet.cnn.layer.SoftmaxLayer'),this.NetworkOutputLayers))
                        % Dissemble dlnetwork to layerGraph to add softmaxLayer
                        % then reconstruct dlnetwork
                        tempRepresentation = layerGraph(this.rlFunction);
                        SoftmaxLayerName = 'RepresentationSoftMax';
                        SoftmaxLayer = softmaxLayer('Name',SoftmaxLayerName);
                        tempRepresentation = addLayers(tempRepresentation,SoftmaxLayer);
                        tempRepresentation = connectLayers(tempRepresentation,this.NetworkOutputNames{1},SoftmaxLayerName);
                        this.rlFunction = dlnetwork(tempRepresentation);
                        this.NetworkOutputNames{1} = SoftmaxLayer.Name;
                        % REVISIT: Update action names, assume discrete
                        % actor representation
                        this.ActionNames{1} = SoftmaxLayer.Name;
                    end
                end
            end
            switch CurrentLoss
                case 'mse'
                    % LossVariable{1} is MSE Target
                    this.LossFcn = @(NetworkOutput, LossVariable) mse(NetworkOutput, LossVariable{1});
                    this.Loss.Param = [];
                case 'cte'
                    % LossVariable{1} is CTE Target
                    this.LossFcn = @(NetworkOutput, LossVariable) rl.loss.cte(NetworkOutput, LossVariable{1}, EntropyLossWeight);
                    this.Loss.Param.EntropyLossWeight = EntropyLossWeight;
                case 'ppo-clipped'
                    this.LossFcn = @(NetworkOutput, LossVariable) rl.loss.ppoClipped(NetworkOutput, LossVariable, EntropyLossWeight, ClipFactor);
                    this.Loss.Param.EntropyLossWeight = EntropyLossWeight;
                    this.Loss.Param.ClipFactor = ClipFactor;
            end
            this.Loss.Name = CurrentLoss;
        end
        
        %% Saveobj
        function s = saveobj(this)
            s.Options = this.Options;
            % dlnetwork will always saved with cpu params
            s.rlFunction = this.rlFunction;
            s.Validated = this.Validated;
            s.ObservationInfo = this.ObservationInfo;
            s.ActionInfo = this.ActionInfo;
            s.ObservationNames = this.ObservationNames;
            s.ActionNames = this.ActionNames;
            s.Loss = this.Loss;
        end
        
        %% Get external network
        function mdl = getModel(this)
            mdl = this.rlFunction;
        end
    end
    
    methods (Static)
        %% Loadobj
        function this = loadobj(s)
            % Recontstruct object
            if isempty(s.ActionInfo)
                this = rlRepresentation(s.rlFunction, s.Options, ...
                    s.ObservationInfo, ...
                    'Observation', s.ObservationNames);
            else
                this = rlRepresentation(s.rlFunction, s.Options, ...
                    s.ObservationInfo, s.ActionInfo, ...
                    'Observation', s.ObservationNames,...
                    'Action', s.ActionNames);
            end
            if ~isempty(s.Loss)
                LossName  = s.Loss.Name;
                LossParam = s.Loss.Param;
                % Set loss
                switch LossName
                    case 'mse'
                        this = setLoss(this,LossName);
                    case 'cte'
                        this = setLoss(this,LossName,'EntropyLossWeight',LossParam.EntropyLossWeight);
                    case 'ppo-clipped'
                        this = setLoss(this,LossName,...
                            'EntropyLossWeight',LossParam.EntropyLossWeight,...
                            'ClipFactor',LossParam.ClipFactor);
                end
            end
        end
    end
    
    %% Implementation of Public Abstract Methods
    methods (Hidden)
        function [gradVal,paramVal] = gradient(this,dy,dx,InputData,varargin)
            % Compute gradients of representation output/loss wrt learnable parameters/inputs
            % This function only supports the following pairs and API
            % dy is 'loss', dx is 'parameter'
            % REVISIT the following when dlnetwork supports MIMO DAG
            % - (dLoss/DParameter,dOutput/dParameter,dOutput/dAction,dOutput/dObservation)
            % - [Y1, ..., YN] = forward(net, X, 'Outputs', layerNames) syntax. 
            %Use case: dOutput/dAction (use name instead of index?)
            
            dx = lower(dx);
            % Map to network inputs
            InputData = mapInputValues(this,InputData);
            % Reshape and cast data to dlarray
            InputData = convertInputDLArray(this, InputData);
            % Cast to gpuArray if GPU training is enabled
            InputData = cellfun(@(x) processData(this,x),InputData,'UniformOutput',false);
            if strcmpi(dy,'output')
                switch dx
                    case 'parameter'
                        % initial gradients only applicable to dOut/dParams
                        if nargin<6
                            initialGradients = [];
                        else
                            initialGradients = varargin{2};
                        end
                        gradVal = gradientOfOutputWrtParams(this,InputData,initialGradients);
                    case 'observation'
                        gradVal = gradientOfOutputWrtInput(this,InputData,this.ObservationIndex);
                    case 'action'
                        if isempty(this.ActionIndex) && ~strcmpi(this.ActionIndex{1},'input')
                            % REVISIT: MSG
                            error('Cannot compute gradient with respect to action since action is not provided as part of the network.')
                        else
                            gradVal = gradientOfOutputWrtInput(this,InputData,this.ActionIndex{2});
                        end
                    otherwise
                        error(message('rl:agent:errInvalidSyntaxGradient'))
                end
            elseif strcmpi(dy,'loss') && strcmpi(dx,'parameter')
                LossVariable = varargin{1};
                gradVal = gradientOfLossWrtParams(this, InputData, LossVariable);
            else
                error(message('rl:agent:errInvalidSyntaxGradient'))
            end
            paramVal = [];
        end
        
        function this = step(this,Gradients,~)
            % Update learnable parameters with gradients with optimizer
            
            % only update params if LR > 0
            if this.Options.LearnRate > 0
                Learnables = getLearnableParameterValues_(this);
                % L2 regularize gradients
                Gradients = rl.internal.dl.regularizeGradient(Gradients, ...
                    Learnables, this.LocalOptions.L2Regularization);
                % Threshold gradients
                Gradients = rl.internal.dl.thresholdGradient(Gradients, ...
                    this.LocalOptions.GradientThreshold, ...
                    this.LocalOptions.GradientThresholdMethod);
                % Parameters update
                % NOTE: ignore LocalLearnRate property in each learnable
                % layers in dlnetwork's Layers until g2017894 is completed
                UpdatedParam = calculateUpdate(this.RLSolver,Learnables,Gradients,this.LocalOptions.InitialLearnRate);
                this = setLearnableParameterValues_(this,UpdatedParam);
            end
        end
        
        function this = fit(this,InputData,LossVariable)
            % Fit learnable parameters on the entire data set
            % Compute dLoss/dParams and update parameters
            
            % only compute gradients and update params if LR > 0
            if this.Options.LearnRate > 0
                Gradients = gradient(this,'loss','parameter',InputData,LossVariable);
                this = step(this,Gradients);
            end
        end
        
        function Representation2 = copy(obj)
            % Copy representation
            Representation2 = obj;
        end
    end
    
    %% Implementation of Protected Abstract Methods
    methods (Access = protected)
        function argStruct = generateEvaluateFunction_(this,argStruct)
            % NOTE: create dummy loss layer and reuse code from rlLayerRepresentation
            
            % generate the evaluate function for the layer representation
            outputstr = argStruct.EvaluateOutputString;
            inputstr  = argStruct.EvaluateInputString;
            policyname = argStruct.PolicyName;
            matfilename = argStruct.MATFileName;
            localfcnstr = argStruct.LocalFunctionString;
            evalfcnname = argStruct.EvaluateFunctionName;
            
            % Extract vars and make sure the DAG can be constructed
            oinfo = this.ObservationInfo;
            ainfo = this.ActionInfo;
            
            if ~isActor(this) || ~isscalar(ainfo) || ~isscalar(oinfo)
                error(message('rl:general:RepLayerCodeGenMIMONotSupported'));
            end
            
            % NOTE: The only place different from rlLayerRep codegen
            % Extract the DAG from dlnetwork
            dag = getDAGNetwork(this); %#ok<NASGU>
            
            % Save the network to disk
            eval(sprintf('%s = dag;',policyname));
            save(matfilename,policyname);
            
            % generate the local function
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
        
        function Value = evaluate_(this, InputData)
            % Perform forward pass, return normal numerical matrix
            % NOTE: return normal numerical matrix since dlarray does not
            % support many function (e.g. cumsum for probability sampling in PG agent)
            % REVISIT with [Y1, ..., YN] = predict(net, X, 'Outputs', layerNames)
            % syntax. Only useful when dlnetwork supports MIMO
            
            % Map to network inputs
            InputData = mapInputValues(this,InputData);
            % Reshape and cast data to dlarray
            InputData = convertInputDLArray(this, InputData);
            % Cast to gpuArray if GPU training is enabled
            InputData = cellfun(@(x) processData(this,x),InputData,'UniformOutput',false);
            % REVISIT: unwrap cell since dlnetwork does not support MIMO
            % Value = extractdata(forward(this.rlFunction, InputData));
            Value = extractdata(predict(this.rlFunction, InputData{1}));
        end
        
        function Value = getLearnableParameterValues_(this)
            % Get dlnetwork learnable params
            Value = this.rlFunction.Learnables;
        end
        
        function this = setLearnableParameterValues_(this,Val)
            % Set dlnetwork learnable params
            this.rlFunction.Learnables = Val;
        end
        
        function validateLearnableParameterSize(this,Value)
            % Bypass this check since dlnetwork already handle learnables
            % parameters datatype and size
        end
        
        function sz = getSize_(this,identifier)
            switch identifier
                case 'parameter'
                    sz = reshape(cellfun(@(x) size(x),getLearnableParameterValues_(this).Value,'UniformOutput',false),1,[]);
                case 'observation'
                    sz = this.ObservationSizes;
                case 'action'
                    sz = this.ActionSizes;
                case 'output'
                    sz = this.OutputSize;
                case 'input'
                    sz = this.InputSize;
            end
        end
        
        function this = preBuildCheckNetwork(this)
            % NOTE: most of the code are borrowed from
            % rlLayerRepresentation. Differences:
            % - Create layerGraph from dlnetwork before the checks
            % - Remove loss layer check since dlnetwork does not have loss
            % - Wait till validateRepresentation (dummy predict pass) to 
            % get network output size since dlnetwork is similar to a 'function'
            
            % do some checks on the network before it is "built"
            
            % get the digraph
            tempLayerGraph = layerGraph(this.rlFunction);
            dg = extractPrivateDirectedGraph(tempLayerGraph);
            layerNodes = dg.Nodes{:,1};
            
            % make sure the graph is DAG
            if ~isdag(dg)
                error(message('rl:agent:errNotDAG'));
            end
            
            % make sure all layers are connected (single bin)
            bins = conncomp(dg,'Type','weak');
            badbins = bins > 1;
            if any(badbins)
                badlayers = arrayfun(@(x) x.Name,layerNodes(badbins),'UniformOutput',false);
                balayersstr = sprintf('"%s", ',badlayers{:});
                balayersstr(end-1:end) = '';
                balayersstr = ['[',balayersstr,']'];
                error(message('rl:agent:errDisconnectedLayers',balayersstr));
            end
            
            % check for any unsupported layers
            if any(arrayfun(@(x) ...
                    isa(x,'nnet.cnn.layer.SequenceInputLayer') || ...
                    isa(x,'nnet.cnn.layer.Image3DInputLayer') || ...
                    isa(x,'nnet.cnn.layer.LSTMLayer') || ...
                    isa(x,'nnet.cnn.layer.BatchNormalizationLayer') || ...
                    isa(x,'nnet.cnn.layer.BiLSTMLayer'),layerNodes))
                error(message('rl:agent:errInvalidLayers'));
            end
            
            % make sure that Normalization is set to None in imageInputLayer
            inputLayer = layerNodes(arrayfun(@(x) isa(x,'nnet.cnn.layer.ImageInputLayer'),layerNodes));
            if any(arrayfun(@(x) ~strcmpi(x.Normalization,'none'),inputLayer))
                error(message('rl:agent:errImageInputLayerNormalizationNone'));
            end
            
            % make sure the IOs are defined
            % get the adjacency matrix
            adj = dg.adjacency';
            % find the output (sink) nodes by finding nodes that are
            % NOT inputs to any other node.
            yidx = find(~any(adj,1));
            for i = 1:numel(yidx)
                yi = yidx(i);
                outLayer = layerNodes(yi);
                if isa(outLayer,'nnet.internal.cnn.layer.SoftmaxLayer') || isa(outLayer,'nnet.cnn.layer.SoftmaxLayer')
                    yidx(i) = yi;
                end
            end
            % find input (src) nodes by finding nodes that have NO incoming
            % connections
            uidx = ~any(adj,2)';
            
            % get the io names
            fcn = @(x) x.Name;
            InputNames = arrayfun(fcn,layerNodes(uidx),'UniformOutput',false);
            OutputNames = arrayfun(fcn,layerNodes(yidx),'UniformOutput',false);
            
            % get the io sizes
            InputSizes = arrayfun(@(x) x.InputSize,layerNodes(uidx),'UniformOutput',false);
            
            % find observation name, index and get its size
            ObservationNames = this.ObservationNames;
            if isempty(intersect(ObservationNames,InputNames,'stable'))
                error(message('rl:agent:errObservationNames'));
            end
            
            % find action name, index and get its size
            ActionNames = this.ActionNames;
            if ~isempty(ActionNames)
                if isempty(intersect(ActionNames,union(InputNames,OutputNames),'stable'))
                    error(message('rl:agent:errActionNames'));
                end
            end
            
            InputNames = reshape(InputNames,1,[]);
            OutputNames = reshape(OutputNames,1,[]);
            
            % find observation name, index and get its size
            ObservationNames = this.ObservationNames;
            [tmp,~,ix] = intersect(ObservationNames,InputNames,'stable');
            ix = reshape(ix,1,[]);
            if isempty(tmp)
                error(message('rl:agent:errObservationNames'));
            end
            this.ObservationIndex = ix;
            this.InputIndex = ix;
            this.ObservationSizes = InputSizes(ix);
            
            % find action name, index and get its size
            ActionNames = this.ActionNames;
            if ~isempty(ActionNames)
                if isempty(intersect(ActionNames,union(InputNames,OutputNames),'stable'))
                    error(message('rl:agent:errActionNames'));
                end
                
                [tmp,~,ix] = intersect(ActionNames,InputNames,'stable');
                ix = reshape(ix,1,[]);
                if ~isempty(tmp)
                    this.ActionIndex = {'input',ix};
                    this.ActionSizes = InputSizes(ix);
                    this.InputIndex = [this.InputIndex this.ActionIndex{2}];
                end
                
                [tmp,~,ix] = intersect(ActionNames,OutputNames,'stable');
                ix = reshape(ix,1,[]);
                if ~isempty(tmp)
                    this.ActionIndex = {'output',ix};
                    % Since action is output, obtain size from validateRep
                end
            end
            
            % validate number of imageInputLayers is equal to observation/action
            numInputs = numel(InputNames);
            numTotalInfo = numel(this.ObservationInfo);
            if ~isActor(this)
                numTotalInfo = numTotalInfo + numel(this.ActionInfo);
            end
            if numInputs ~= numTotalInfo
                error(message('rl:agent:errIncompatibleNumberOfImageInputLayer',numInputs,numTotalInfo));
            end
            this.InputSize = InputSizes;
            this.NetworkInputNames  = InputNames;
            this.NetworkOutputNames = OutputNames;
            this.NetworkOutputLayers = layerNodes(yidx);
        end
        
        function this = updateOptions(this,opt)
            % convert representation option to local options and update
            % solver
            this.LocalOptions = rl.util.rlLayerRepresentation.createOptimizerOptions(opt);
            switch class(this.LocalOptions)
                case "nnet.cnn.TrainingOptionsADAM"
                    this.RLSolver = rl.internal.dl.rlSolverADAM(this.LocalOptions);
                case "nnet.cnn.TrainingOptionsSGDM"
                    this.RLSolver = rl.internal.dl.rlSolverSGDM(this.LocalOptions);
                case "nnet.cnn.TrainingOptionsRMSProp"
                    this.RLSolver = rl.internal.dl.rlSolverRMSProp(this.LocalOptions);
            end
        end
    end
    
    %% Compute Gradient methods
    methods (Access = private)
        % REVISIT currently dlnetwork does not support MIMO DAG
        % Only support dLoss/dParams for PG,AC,PPO due to MIMO limitation
        function [GradVal, Loss] = gradientOfLossWrtParams(this, InputData, LossVariable)
            [GradVal, Loss] = dlfeval(@dLoss_dParams, this, this.rlFunction, InputData, LossVariable);
        end
        
        function [GradVal, Loss] = dLoss_dParams(this, DLNet, InputData, LossVariable)
            % dLoss/dParams
            % REVISIT: unwrap InputData since dlnetwork does not support MIMO yet
            % NetworkOutput = forward(DLNet, InputData);
            NetworkOutput = forward(DLNet, InputData{1});
            Loss = feval(this.LossFcn, NetworkOutput, LossVariable);
            GradVal = dlgradient(Loss, DLNet.Learnables);
        end
        
        function GradVal = gradientOfOutputWrtInput(this, InputData, InputIdx)
            FcnHandle = @(InputData) dOutput_dInput(this, InputData, InputIdx);
            GradVal = dlfeval(FcnHandle, InputData);
        end
        
        function GradVal = dOutput_dInput(this, InputData, InputIdx)
            % dOut/dIn, REVISIT: How to handle multi action channels
            NetworkOutput = forward(this.Model, InputData, getLearnableParameterValues_(this));
            GradVal = dlgradient(NetworkOutput, InputData{InputIdx});
        end
        
        function GradVal = gradientOfOutputWrtParams(this, InputData, InitialGradient)
            FcnHandle = @(Params) dOut_dParams(this, InputData, Params, InitialGradient);
            GradVal = dlfeval(FcnHandle, getLearnableParameterValues_(this));
        end
        
        function GradVal = dOut_dParams(this, InputData, Params, InitialGradient)
            % dOut/dParams
            if isempty(InitialGradient)
                Output = forward(this.Model, InputData, Params);
            else
                % InitialGradient is similar to dLoss/dOutput
                Output = sum(forward(this.Model, InputData, Params) .* InitialGradient,'all');
            end
            GradVal = dlgradient(Output, Params);
        end
    end
    
    methods (Hidden)
        %% Parameters target update (overload)
        function this = updateParameter(this,Representation,smoothingFactor)
            % REVISIT: this method should be abstract on the abstract rep class
            % Update learnable parameters of one representation with others
            
            if nargin<3
                smoothingFactor = 1;
            end
            
            NewParameter = getLearnableParameterValues_(Representation);
            CurrentParameter = getLearnableParameterValues_(this);
            CurrentParameter = ...
                dlupdate(@(currentParams,newParams) ...
                    currentParams.*(1-smoothingFactor) + newParams.*smoothingFactor,...
                    CurrentParameter, NewParameter);
            
            this = setLearnableParameterValues_(this,CurrentParameter);
        end
        
        % REVISIT: Remove once we completely deprecate rlModelRep
        function Params = getParameters(this)
            Params = getLearnableParameterValues_(this);
        end
        
        %% Validate with dummy predict pass
        function this = validateRepresentation(this)
            % NOTE: borrow input/output index mixing from rlLayerRep
            % Obtain network output size after dummy predict pass since
            % cannot query output size from dlnetwork API
            
            % Number of observationInfo, actionInfo have to match number of names
            if numel(this.ObservationInfo)~=numel(this.ObservationNames)
                error(message('rl:agent:errIncompatibleObservation1'));
            end
            if numel(this.ActionInfo)~=numel(this.ActionNames)
                error(message('rl:agent:errIncompatibleAction1'));
            end
            
            validateRepresentationObservationSizes(this,this.ObservationSizes);
            % if action at output, validate after forward pass
            if ~isempty(this.ActionIndex) && strcmpi(this.ActionIndex{1},'input')
                validateRepresentationActionSizes(this,this.ActionSizes);
            end
            
            % usample from data spec and validate network by predict pass
            try
                DummyInput = usampleInputs(this);
                DummyInput = convertInputDLArray(this, DummyInput);
                % REVISIT: unwrap cell since dlnetwork does not support MIMO
                % ForwardOutput = extractdata(predict(this.rlFunction, DummyInput));
                ForwardOutput = extractdata(predict(this.rlFunction, DummyInput{1}));
                % REVISIT: currently only allow single output channel 
                % (single output action channel or critic)
                this.OutputSize{1} = size(ForwardOutput);
                if isActor(this)
                    this.ActionSizes = this.OutputSize;
                    validateRepresentationActionSizes(this,this.ActionSizes);
                end
                
                this.Validated = true;
            catch ex
                this.Validated = false;
                me = MException(message('rl:agent:errRepresentationEvaluation'));
                throw(addCause(me,ex))
            end
        end
        
        %% Construct DAGNetwork from dlnetwork
        function dag = getDAGNetwork(this)
            % Create layerGraph from dlnetwork
            lg = layerGraph(this.rlFunction);
            % NOTE: Create dummy loss layer to support codegen 19b
            LossName = 'RepresentationLoss';
            switch this.Loss.Name
                case 'mse'
                    LossLayer = regressionLayer('Name', LossName);
                case 'cte'
                    if this.Loss.Param.EntropyLossWeight
                        LossLayer = rl.layer.internal.ClassificationLayerWithEntropy(...
                            'Name',LossName,'EntropyLossWeight',this.Loss.Param.EntropyLossWeight);
                    else
                        LossLayer = classificationLayer('Name', 'RepresentationLoss');
                    end
                case 'ppo-clipped'
                    LossLayer = rl.layer.internal.ClipPPOLossLayer(...
                        'Name',LossName,...
                        'EntropyLossWeight',this.Loss.Param.EntropyLossWeight,...
                        'ClipFactor',this.Loss.Param.ClipFactor);
            end
            if ismember(this.Loss.Name, {'cte','ppo-clipped'})
                % Store classification output meta data (class name)
                ActionValues = reshape(this.ActionInfo.Elements,[],1);
                if iscell(ActionValues)
                    ActionValues = cellfun(@mat2str,ActionValues,'UniformOutput',false);
                end
                LossLayer.Classes = categorical(ActionValues, ActionValues);
            end
            % Attach dummy loss layer to layerGraph
            lg = addLayers(lg, LossLayer);
            lg = connectLayers(lg, this.NetworkOutputNames{1}, LossName);
            % Construct DAG network from layerGraph
            dag = assembleNetwork(lg);
        end
        
        %% Get loss
        function Loss = getLoss(this)
            Loss = this.Loss;
        end
    end
    
    %% Data manipulation
    methods (Access = private)
        function Data = convertInputDLArray(this, Data)
            % Reshape and cast numeric input data to dlarray with correct
            % dimension label.
            % REVISIT: handle dlarray input, when customer write custom
            % train loop
            
            if ~iscell(Data)
                Data = {Data};
            end
            
            % Find batch size by dividing the size of data from 1 channel
            % to input size from that channel.
            % NOTE: Cheap computation since numel() uses metadata of matrix 
            BatchSize = numel(Data{1}) / prod(this.InputSize{1});
            
            % REVISIT: for 19b, only support imageInputLayer (SSCB)
            for ix = 1:numel(Data)
                if BatchSize > 1
                    switch ndims(Data{ix})
                        case 4
                            Data{ix} = dlarray(Data{ix}, 'SSCB');
                        case {3, 2}
                            Data{ix} = dlarray(Data{ix}, [repmat('S',1,numel(this.InputSize{ix})) 'B']);
                        otherwise
                            error(message('rl:agent:errDLNetworkInputDimNotSupport'))
                    end
                else
                    switch ndims(Data{ix})
                        case {3, 2, 1}
                            Data{ix} = dlarray(Data{ix}, 'SSC');
                        otherwise
                            error(message('rl:agent:errDLNetworkInputDimNotSupport'))
                    end
                end
            end
        end
        
        function Data = processData(this,Data)
            % Cast to gpuArray if GPU training is enabled
            % Modify from rlLayerRepresentation
            
            if useGPU(this)
                Data = gpuArray(Data);
            else
                Data = single(Data);
            end
        end
        
        function flag = useGPU(this)
            flag = ismember(this.LocalOptions.ExecutionEnvironment, {'gpu'} );
        end
    end
end