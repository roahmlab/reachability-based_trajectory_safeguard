%rlAbstractSamplingStrategy Defines the interface of stochastic policy 
% sampling strategy.
%
% rlAbstractSamplingStrategy public methods
%    sample   - Return a sample or batch samples from given probability
%               distribution defined by ProbabilityFcn and input 
%               ProbabilityParameters.
%    evaluate - Evaluate probability function, used by autodiff backward.

% Copyright 2019 The MathWorks, Inc.

classdef rlAbstractSamplingStrategy
    
    properties (Access = protected)
        SamplingSpec
    end
    
    methods
        function this = rlAbstractSamplingStrategy(SamplingSpec)
            % Constructor
            
            validateattributes(SamplingSpec, {'rl.util.RLDataSpec'}, {'vector', 'nonempty'}, 'SamplingSpec');
            this.SamplingSpec = SamplingSpec;
        end
        function Out = sample(this, varargin)
            %SAMPLE return a sample or batch samples from probability
            % distribution defined by ProbabilityFcn and input 
            % ProbabilityParameters.
            
            Out = sampleImpl(this, varargin{:});
        end
        
        function Out = evaluate(this, varargin)
            %EVALUATE evaluate the policy density function. e.g. pi(x|mu,sigma)
            % ActionProbabilityParam (e.g. mean, standard deviation) at Action (x)
            %   Policy = ActionProbabilityParam (action probabilities) for 
            %   discrete action
            
             Out = distributionFcn(this, varargin{:});
        end
        
        function Model = setModelOutputType(this, Model)
            %Model = SETMODELOUTPUTTYPE(SamplingStrategy, Model) set and
            % validate the model output type to be compatible with the
            % sampling strategy.
            % E.g. discrete sampling strategy will add softmax to the end
            % of the model output (output probability)
            
            Model = setModelOutputTypeImpl(this, Model);
        end
    end
    
    methods (Hidden)
        % REVISIT: make hidden since this is internal operation
        function validateModelOutputCompatibility(this, ModelOutputSize)
            %VALIDATEMODELOUTPUTCOMPATIBILITY checks if model output size
            % is compatible with the sampling strategy. e.g. Continuous
            % Gaussian strategy requires model output size of 2*nAction
            
            % implement of each class
            % if sampling strategy == none, check num output = num action
            % if sampling strategy == gaussian, check output size (mean + standard deviation)
            % if sampling strategy == discrete, check output size (action probability)
        end
    end
    
    methods (Access = protected)
        out = sampleImpl(this, varargin);
        %SAMPLEIMPL return a sample or batch samples from probability
        % distribution defined by ProbabilityFcn and input
        % ProbabilityParameters.
        
        out = distributionFcn(this, varargin);
        %EVALUATEIMPL evaluate the policy density function. e.g. pi(x|mu,sigma)
        % ActionProbabilityParam (e.g. mean, standard deviation) at Action (x)
        %   Policy = ActionProbabilityParam (action probabilities) for
        %   discrete action
        
        Model = setModelOutputTypeImpl(this, Model);
        %Model = SETMODELOUTPUTTYPEIMPL(SamplingStrategy, Model) set and
        % validate the model output type to be compatible with the
        % sampling strategy.
        % E.g. discrete sampling strategy will add softmax to the end
        % of the model output (output probability)
    end
end