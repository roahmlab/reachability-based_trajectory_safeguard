classdef rlModelRepresentation < rl.util.rlAbstractRepresentation
% rlModelRepresentation: RL representation for RL agents using a neural network
% from rl.ad library
%
% To create an rlModelRepresention object, use the rlRepresentation
% function.
%
%   rep = rlRepresentation(ADMODEL,OINFO,'Observation',ONAMES) creates a
%   representation with default options using ADMODEL, which is a neural
%   network from the rl.ad library. Specify the network input layer names,
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
        DimensionsToSqueeze
        ADVersion
        
        % model io names determined before compile()
        ModelInputNames
        ModelOutputNames
    end
    
    methods
        %% Constructor
        function this = rlModelRepresentation(rlFunction,varargin)
            this = this@rl.util.rlAbstractRepresentation(rlFunction,varargin{:});
            
            if isa(rlFunction,'rl.internal.ad.Model')
                this.ADVersion = 1;
            else
                this.ADVersion = 2;
            end
            this.adModel = rlFunction;
            this.LocalOptions = rl.util.rlModelRepresentation.createOptimizerOptions(this.Options,this.ADVersion);
            setOptimizer(this.adModel,this.LocalOptions);
            
            this.ModelInputNames  = this.adModel.getInputs;
            this.ModelOutputNames = this.adModel.getOutputs;
                        
            if this.adModel.IsBuilt
                this = validateRepresentation(this);                
            end
        end
        
        %% Set loss
        function this = setLoss(this,AutomaticLoss,varargin)
            AutoLoss = validatestring(AutomaticLoss,{'mse','cte'},'','AutomaticLoss');
            
            % Input parser
            parser = inputParser;
            addParameter(parser,'EntropyLossWeight',0);
            parse(parser,varargin{:});
            r = parser.Results;
            EntropyLossWeight = r.EntropyLossWeight;
            validateattributes(EntropyLossWeight,{'numeric'},{'scalar','nonnegative','finite','real','<=',1},'','EntropyLossWeight');
            
            % Validates whether network is compiled and compiles network if necessary
            rlFunction = this.adModel;
            
            if ~rlFunction.AutomaticCompilationRequired && ~rlFunction.IsBuilt
                error(message('rl:agent:errModelOptimizerNeeded'));
            end
            if rlFunction.AutomaticCompilationRequired
                if strcmpi(AutoLoss,'cte')
                    % append cte loss
                    if ~isa(this.ActionInfo,'rl.util.rlFiniteSetSpec')
                        error(message('rl:agent:errCTEDiscreteOutput'));
                    end
                    % check and append softmax
                    CurrentInputNode = rlFunction.getInputs;
                    % TODO: v1 supports only single loss
                    CurrentOutputNode = rlFunction.getOutputs{1};
                    
                    if ~(isa(CurrentOutputNode,'rl.internal.ad.ops.Softmax') || isa(CurrentOutputNode,'rl.ad.ops.Softmax'))
                        % reconstruct adModel: wrap output with softmax
                        if this.ADVersion == 1
                            rlFunction = rl.internal.ad.Model(CurrentInputNode,softmax(CurrentOutputNode));
                        else
                            rlFunction = rl.ad.Model(CurrentInputNode,softmax(CurrentOutputNode));
                        end
                        setOptimizer(rlFunction,this.LocalOptions);
                    end
                end
                compile(rlFunction,AutoLoss,rlFunction.LocalOptimizer,'EntropyLossWeight',EntropyLossWeight);
            else
                error(message('rl:agent:errLossDefined'));
            end
            this.adModel = rlFunction;
            this = validateRepresentation(this);
        end        
        
        %% Get model
        function mdl = getModel(this)
            % REVISIT bypass for now
        end
    end
    methods (Hidden)
        function Parameters = getParameters(this)
            Parameters = getLearnableParameters(this.adModel);
        end

        %% CHECK: Compute gradients of representation output wrt learnable parameters or inputs
        function [grad,parameter] = gradient(this,dy,dx,inputValues,varargin)
            % This function only supports the following pairs and API
            %     dy       dx
            %   'loss'   'param'  : gradVal = gradient(rep,'loss','parameter',InputData,Target);
            %  'output'  'action' : gradVal = gradient(rep,'output','action',InputData);
            %  'output'  'param'  : gradVal = gradient(rep,'output','action',InputData,'initialGradients',initGradVal);            
            
            if nargin<6
                initialGrad = false;
            else
                initialGrad = true;
            end
            
            % TODO: What happens when actions, observations, targets, outputs are MIMO?
            % TODO: What happens there is multiple output
            
            % get output and input nodes            
            dy = lower(dy);
            dx = lower(dx);
            switch lower(dy)
                case 'output'
                    dOutput = getOutputs(this.adModel,1);
                case 'loss'
                    dOutput = getLoss(this.adModel);
                otherwise
                    error(message('rl:agent:errUnknownIdentifier'));
            end
            switch lower(dx)
                case 'parameter'
                    dInput = [];
                case 'observation'
                    dInput = getInputs(this.adModel,this.ObservationIndex);
                case 'action'
                    dInput = getInputs(this.adModel,this.ActionIndex{2});
                otherwise
                    error(message('rl:agent:errUnknownIdentifier'));
            end
            
            % set input values
            Data = mapInputValues(this,inputValues);
            Data = cellfun(@(x) processData(this,x), Data, 'UniformOutput', false);
            Data = convertInputData(this,Data);
            Inputs = getInputs(this.adModel);            

            for ct=1:length(Inputs)
                Inputs{ct}.Value = Data{ct};
            end            
            
            % compute gradients
            if strcmpi(dx,'observation') || strcmpi(dx,'action')
                grad = gradient(dOutput,dInput);
                parameter = [];
            elseif strcmpi(dx,'parameter')
                optimizer = getOptimizer(this.adModel);
                if nargin == 4
                    [grad,parameter] = computeGradients(optimizer,...
                        dOutput,dInput);
                else
                if initialGrad
                    [grad,parameter] = computeGradients(optimizer,...
                        dOutput,dInput,varargin{1},varargin{2});
                else
                    % set input values
                    Targets = getTargets(this.adModel);
                    for ct=1:length(Targets)
                         Target = processData(this,varargin{1}{1});
                         Targets{ct}.Value = convertTargetData(this,Target);
                    end                    
                    [grad,parameter] = computeGradients(optimizer,dOutput);
                end
                end
            end
        end
        
        %% Step learnable parameters with input and target data        
        function this = fit(this,Data,Target)
            % TODO: What happens when data, and targets are MIMO?
            Data = mapInputValues(this,Data);
            Data = cellfun(@(x) processData(this,x),Data,'UniformOutput',false);
            Data = convertInputData(this,Data);
            % TODO: Determine BatchSize reliably
            sz = size(Target{1});
            BatchSize = sz(end);     
            Target = cellfun(@(x) processData(this,x),Target,'UniformOutput',false);
            Target = convertTargetData(this,Target);
            if this.ADVersion == 1
                Batch = rl.internal.ad.datasource.xyDataSource(Data,Target,BatchSize);
            else
                Batch = rl.ad.datasource.xyDataSource(Data,Target,BatchSize);
            end
            fit(this.adModel,Batch,1,isGPU(this),'Verbose',0);
        end
        
        %% Step learnable parameters with gradient information or loss
        function this = step(this,gradient,parameter)
            optimizer = getOptimizer(this.adModel);
            applyGradients(optimizer,parameter,gradient);
        end                
        
        %% Copy representation
        function Representation2 = copy(Representation1)
            % Representation is value object
            rlFunction = copy(Representation1.adModel);
            % Create new one
            Representation2 = rlRepresentation(rlFunction,Representation1.Options, ...
                                     'Observation',Representation1.ObservationNames,Representation1.ObservationInfo, ...
                                     'Action',Representation1.ActionNames,Representation1.ActionInfo);
        end
    end
    
    methods (Access = protected)
       
        function Value = evaluate_(this,Data)
            % Map to network inputs
            Data = mapInputValues(this,Data);
            Data = cellfun(@(x) processData(this,x), Data, 'UniformOutput', false);            
            Data = convertInputData(this,Data);
            Value = predict(this.adModel,Data,isGPU(this));
            Value = convertOutputData(this,Value);
        end
        function Value = getLearnableParameterValues_(this)
            params = getParameters(this);
            for ct=1:numel(params)
                Value{1,ct} = params{ct}.Value;                
            end
        end
        function this = setLearnableParameterValues_(this,Value)
            params = getParameters(this);
            for ct=1:length(params)
                params{ct}.Value = processData(this,Value{ct});
            end            
        end
        function sz = getSize_(this,identifier)
            switch identifier
                case 'observation'
                    Inputs = getInputs(this.adModel);
                    ix = this.ObservationIndex;
                    for ct=1:length(ix)
                        sz{1,ct} = Inputs{ix(ct)}.Size;
                    end
                case 'action'
                    sz = {};
                    if ~isempty(this.ActionNames)
                        ix = this.ActionIndex{2};
                        if strcmpi(this.ActionIndex{1},'input')
                            Nodes = getInputs(this.adModel);
                        else
                            Nodes = getOutputs(this.adModel);                            
                        end
                        for ct=1:length(ix)
                            sz{1,ct} = Nodes{ix(ct)}.Size;
                        end                        
                    end
                case 'output'                    
                    Outputs = getOutputs(this.adModel);
                    for ct=1:numel(Outputs)
                        sz{1,ct} = Outputs{ct}.Size;                        
                    end
                case 'input'
                    Inputs = getInputs(this.adModel);
                    for ct=1:numel(Inputs)
                        sz{1,ct} = Inputs{ct}.Size;
                    end
                case 'parameter'
                    params = getParameters(this);
                    sz = cellfun(@(x) x.Size,params,'UniformOutput',false);
            end
        end
        function this = validateRepresentation(this)
            Inputs  = this.ModelInputNames;
            Outputs = this.ModelOutputNames;
            
            InputNames = cellfun(@(x) x.Name,Inputs,'UniformOutput',false);
            OutputNames = cellfun(@(x) x.Name,Outputs,'UniformOutput',false);
            
            %% find observation name, index and get its size
            ObservationNames = this.ObservationNames;
            [tmp,ix] = intersect(ObservationNames,InputNames,'stable');
            ix = reshape(ix,1,[]);
            if isempty(tmp)
                error(message('rl:agent:errObservationNames'));
            end
            this.ObservationIndex = ix;
            ObservationSizes = cellfun(@(x) x.Size,Inputs(ix),'UniformOutput',false);
            this.InputIndex = ix(:)';
            
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
                    ActionSizes = cellfun(@(x) x.Size,Inputs(ix),'UniformOutput',false);
                    this.InputIndex = [this.InputIndex this.ActionIndex{2}];
                end
                
                [tmp,~,ix] = intersect(ActionNames,OutputNames,'stable');
                ix = reshape(ix,1,[]);
                if ~isempty(tmp)
                    this.ActionIndex = {'output',ix};
                    ActionSizes = cellfun(@(x) x.Size,Outputs(ix),'UniformOutput',false);
                end                                
            end
                
            validateRepresentationSizes(this,ObservationSizes,ActionSizes);
            
            % define dimensions to squeeze
            InputSizes = cellfun(@(x) x.Size,Inputs,'UniformOutput',false);
            this.DimensionsToSqueeze = false(size(InputSizes));
            for ct=1:numel(InputSizes)
                sz = InputSizes{ct};
                if length(sz)==2
                    if sz(2)==1
                       this.DimensionsToSqueeze(ct) = true; 
                    end
                else
                    % size == 3
                    if sz(2)==1 && sz(3)==1
                        this.DimensionsToSqueeze(ct) = true;
                    end
                end              
            end
            
            this.Validated = true;
        end        
        function validateLearnableParameterSize(this,Value)
            if ~isa(Value,'cell')
                error(message('rl:agent:errLearnableParameters'));
            end
            
            params = getLearnableParameterValues(this);
            paramSize = cellfun(@(x) size(x),params,'UniformOutput',false);
            % REVIST for some reason the size propery of the parameters
            % goes stale ([]) on workers resulting in parallel failures
            
            % params = getParameters(this);
            % paramSize = cellfun(@(x) x.Size,params,'UniformOutput',false);
            valSize = cellfun(@(x) size(x),Value,'UniformOutput',false);
            match = cellfun(@(x,y) isequal(x,y),paramSize,valSize);
            if ~all(match)
                idx = find(~match,1);
                error(message('rl:agent:errParameterMismatch',mat2str(paramSize{idx}),mat2str(valSize{idx})));
            end            
        end        
        function Data = convertInputData(this,Data)
            dt = Data(this.DimensionsToSqueeze);
            for ct=1:numel(dt)
                sz = size(dt{ct},1);
                dt{ct} = reshape(dt{ct},sz,[]);                
            end
            Data(this.DimensionsToSqueeze) = dt;
        end
        function Data = convertOutputData(this,Data)
            sz = size(Data);
            sz = [sz ones(1,4-length(sz))];
            Data = reshape(Data,[sz(1) 1  1 sz(2)]);
            
            if isGPU(this)
                Data = gather(Data);
            end
        end        
        function Data = convertTargetData(this,Data)
            if iscell(Data)
                for ct=1:numel(Data)
                    sz = size(Data{ct});
                    Data{ct} = reshape(Data{ct},[sz(1) sz(end) 1 1]);
                end
            else
                sz = size(Data);
                Data = reshape(Data,[sz(1) sz(end) 1 1]);                
            end
        end
        function this = updateOptions(this,opt)
            % convert representation option to local options
            this.LocalOptions = rl.util.rlModelRepresentation.createOptimizerOptions(opt,this.ADVersion);
            if this.Validated                
                setOptimizer(this.adModel,this.LocalOptions);                
            end
        end
    end
    methods (Access = private)
        function Data = processData(this,Data)
            if isGPU(this)
                Data = gpuArray(single(Data));
            else
                Data = single(Data);
            end
        end
        function flag = isGPU(this)
            flag = strcmpi(this.Options.UseDevice,'GPU');
        end
    end
    methods (Static,Hidden)
        function newopt = createOptimizerOptions(opt,ADVersion)
            % supported optimizers
            optimizer = validatestring(opt.Optimizer,{'adam','sgdm','rmsprop'},'','Optimizer');
            if ADVersion == 1
                switch optimizer
                    case 'adam'
                        newopt = rl.internal.ad.optimizer.Adam('LearnRate',opt.LearnRate, ...
                            'Beta1',opt.OptimizerParameters.GradientDecayFactor, ...
                            'Beta2',opt.OptimizerParameters.SquaredGradientDecayFactor, ...
                            'Epsilon',opt.OptimizerParameters.Epsilon, ...
                            'GradientThresholdMethod',opt.GradientThresholdMethod, ...
                            'GradientThreshold',opt.GradientThreshold, ...
                            'L2Regularization',opt.L2RegularizationFactor);
                    case 'sgdm'
                        newopt = rl.internal.ad.optimizer.Sgdm('LearnRate',opt.LearnRate, ...
                            'Momentum',opt.OptimizerParameters.Momentum, ...
                            'GradientThresholdMethod',opt.GradientThresholdMethod, ...
                            'GradientThreshold',opt.GradientThreshold, ...
                            'L2Regularization',opt.L2RegularizationFactor);
                    case 'rmsprop'
                        newopt = rl.internal.ad.optimizer.Rmsprop('LearnRate',opt.LearnRate, ...
                            'Beta2',opt.OptimizerParameters.SquaredGradientDecayFactor, ...
                            'Epsilon',opt.OptimizerParameters.Epsilon, ...
                            'GradientThresholdMethod',opt.GradientThresholdMethod, ...
                            'GradientThreshold',opt.GradientThreshold, ...
                            'L2Regularization',opt.L2RegularizationFactor);
                end
            else
                switch optimizer
                    case 'adam'
                        newopt = rl.ad.optimizer.Adam('LearnRate',opt.LearnRate, ...
                            'Beta1',opt.OptimizerParameters.GradientDecayFactor, ...
                            'Beta2',opt.OptimizerParameters.SquaredGradientDecayFactor, ...
                            'Epsilon',opt.OptimizerParameters.Epsilon, ...
                            'GradientThresholdMethod',opt.GradientThresholdMethod, ...
                            'GradientThreshold',opt.GradientThreshold, ...
                            'L2Regularization',opt.L2RegularizationFactor);
                    case 'sgdm'
                        newopt = rl.ad.optimizer.Sgdm('LearnRate',opt.LearnRate, ...
                            'Momentum',opt.OptimizerParameters.Momentum, ...
                            'GradientThresholdMethod',opt.GradientThresholdMethod, ...
                            'GradientThreshold',opt.GradientThreshold, ...
                            'L2Regularization',opt.L2RegularizationFactor);
                    case 'rmsprop'
                        newopt = rl.ad.optimizer.Rmsprop('LearnRate',opt.LearnRate, ...
                            'Beta2',opt.OptimizerParameters.SquaredGradientDecayFactor, ...
                            'Epsilon',opt.OptimizerParameters.Epsilon, ...
                            'GradientThresholdMethod',opt.GradientThresholdMethod, ...
                            'GradientThreshold',opt.GradientThreshold, ...
                            'L2Regularization',opt.L2RegularizationFactor);
                end
            end
        end
        function validateRepresentationType(rep)
            if ~isa(rep,'rl.util.rlModelRepresentation')
                error(message('rl:agent:errInvalidRepresentation'));
            end
        end         
    end
end