classdef rlContinuousGaussianStrategy < rl.representation.sampling.rlAbstractContinuousStrategy
    % Continuous uniform sampling object
    
    % Copyright 2019 The MathWorks, Inc.
    
    properties (Access = private)
        % Number of action, calculated from data spec.
        % REVISIT: should this be on abstract continuous strategy
        NumAction
        NumActionChannel
        ActionIdxPerChannel
        
        MeanIdx
        StandardDeviationIdx
    end
    
    methods
        function this = rlContinuousGaussianStrategy(SamplingSpec)
            
            this = this@rl.representation.sampling.rlAbstractContinuousStrategy(SamplingSpec);
            this.NumActionChannel = numel(this.SamplingSpec);
            NumActionPerChannel = cellfun(@prod, {this.SamplingSpec.Dimension});
            this.NumAction = sum(NumActionPerChannel);
            this.MeanIdx = 1 : this.NumAction;
            this.StandardDeviationIdx = (this.NumAction+1) : (this.NumAction*2);
            
            % Map to action info
            this.ActionIdxPerChannel = cell(this.NumActionChannel,1);
            LastIdx = 0;
            AllActionIdx = 1:this.NumAction;
            for ct = 1:this.NumActionChannel
                NewIdx = LastIdx + NumActionPerChannel(ct);
                this.ActionIdxPerChannel{ct} = AllActionIdx(LastIdx+1:NewIdx);
                LastIdx = NewIdx;
            end
        end
    end
    
    methods (Access = protected)
        function Action = sampleImpl(this, ProbabilityParam, BatchSize, SequenceLength)
            % Return random sample from continuous normal distribution
            % with specific mean and variance. Use in getAction(rlActorRepresentation).
            
            % NOTE: only support is single channel output
            ProbabilityParam = squeeze(ProbabilityParam{1});
            DimsToSlice = 1;
            Mean = rl.internal.dataTransformation.generalSubref(ProbabilityParam,this.MeanIdx,DimsToSlice);
            StandardDeviation = rl.internal.dataTransformation.generalSubref(ProbabilityParam,this.StandardDeviationIdx,DimsToSlice);
            
            if any(StandardDeviation < 0,'all')
                error(message('rl:agent:errNegativeStandardDeviation'));
            end
            
            Out = Mean + StandardDeviation .* randn(size(Mean),'like',Mean);
            Action = cell(1,numel(this.ActionIdxPerChannel));
            % Map to action info
            for ct = 1:this.NumActionChannel
                Action{ct} = rl.internal.dataTransformation.generalSubref(Out,this.ActionIdxPerChannel{ct},DimsToSlice);
            end
        end
        
        function Out = distributionFcn(this, varargin)
            % Normal probability density function. Use for backpropagation.
            
            ProbabilityParam = squeeze(varargin{1}); % mean and variance
            ValueAtSample = varargin{2};             % sampled action

            % after squeeze(), always slice on 1st dimension
            DimsToSlice = 1;
            Mean = rl.internal.dataTransformation.generalSubref(ProbabilityParam,this.MeanIdx,DimsToSlice);
            StandardDeviation = rl.internal.dataTransformation.generalSubref(ProbabilityParam,this.StandardDeviationIdx,DimsToSlice);
            ValueAtSample = reshape(ValueAtSample,size(Mean));
            
            if any(StandardDeviation < 0,'all')
                error(message('rl:agent:errNegativeStandardDeviation'));
            end
            
            % protect divide by zeros
            StandardDeviation = rl.internal.dataTransformation.boundAwayFromZero(StandardDeviation);
            Out = rl.internal.dataTransformation.boundAwayFromZero(exp(-0.5 * ((ValueAtSample - Mean)./StandardDeviation).^2)) ./ (sqrt(2*pi) .* StandardDeviation);
        end
        
        function Model = setModelOutputTypeImpl(this, Model)
            %Model = SETMODELOUTPUTTYPEIMPL(SamplingStrategy, Model) set/
            % validate the model output type to be compatible with the
            % sampling strategy.
            % gaussian sampling strategy will add tanh and scale the output
            % for the mean path according to sampling spec, add softplus to 
            % the variance path (ensure positive)
            
            [ModelOutputTransform.Scale, ModelOutputTransform.Bias] = getTanhScaleBias(this.SamplingSpec);
            PostProcessFcn = @(ModelOutput) processModelOutput(this,ModelOutput,ModelOutputTransform);
            Model = setOutputType(Model,'gaussian',PostProcessFcn);
        end
    end
    
    methods (Access = private)
        function ModelOutput = processModelOutput(this, ModelOutput, ModelOutputTransform)
            
            DimsToSlice = 1;
            Mean = rl.internal.dataTransformation.generalSubref(ModelOutput,this.MeanIdx,DimsToSlice);
            StandardDeviation = rl.internal.dataTransformation.generalSubref(ModelOutput,this.StandardDeviationIdx,DimsToSlice);
            
            % bound and shift the Mean to value of interest with
            % tanh + transform
            ModelOutput(this.MeanIdx,:,:) = tanh(Mean) * ModelOutputTransform.Scale{1} + ModelOutputTransform.Bias{1};
            
            % apply softplus to StandardDeviation
            ModelOutput(this.StandardDeviationIdx,:,:) = rl.layer.SoftplusLayer.evaluate(StandardDeviation);
        end
    end
end