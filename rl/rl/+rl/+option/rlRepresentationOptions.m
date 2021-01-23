classdef rlRepresentationOptions
    % rlRepresentationOptions: Creates options for agent representations.
    
    % Copyright 2017-2018 The MathWorks Inc.
    
    properties
        LearnRate
        GradientThreshold
        GradientThresholdMethod        
        L2RegularizationFactor  
        UseDevice
    end
    properties (Dependent)
        Optimizer
        OptimizerParameters
    end
    properties (Access = private)
        Optimizer_
        OptimizerParameters_
    end
	properties (Hidden)
        AllowMultipleActionChannels (1,1) logical = false
    end
    methods
        function obj = rlRepresentationOptions(varargin)
            parser = inputParser;
            addParameter(parser,'LearnRate',1e-2);
            addParameter(parser,'Optimizer',"adam");            
            addParameter(parser,'GradientThreshold',Inf);
            addParameter(parser,'GradientThresholdMethod',"l2norm");
            addParameter(parser,'L2RegularizationFactor',1e-4);
            addParameter(parser,'UseDevice',"cpu");
            
            parse(parser,varargin{:});
            
            obj.LearnRate = parser.Results.LearnRate;
            obj.Optimizer = parser.Results.Optimizer;
            obj.GradientThreshold = parser.Results.GradientThreshold;
            obj.GradientThresholdMethod = parser.Results.GradientThresholdMethod;
            obj.L2RegularizationFactor = parser.Results.L2RegularizationFactor;
            obj.UseDevice = parser.Results.UseDevice;
        end
        function obj = set.LearnRate(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','nonnegative','finite'},'','LearnRate');
            obj.LearnRate = Value;
        end
        function obj = set.Optimizer(obj,Value)
            validstr = validatestring(Value,{'adam','sgdm','rmsprop'},'','Optimizer');
            obj.Optimizer_ = string(validstr);
            if isempty(obj.OptimizerParameters_)
                obj.OptimizerParameters_ = rl.option.OptimizerParameters(obj.Optimizer);
            else
                obj.OptimizerParameters_ = setOptimizer(obj.OptimizerParameters,obj.Optimizer);
            end
        end
        function obj = set.OptimizerParameters(obj,Value)
            validateattributes(Value,{'rl.option.OptimizerParameters'},{'scalar'},'','OptimizerParameters')
            obj.OptimizerParameters_ = Value;
            NewOptimizer = getOptimizer(Value);
            if ~isequal(obj.Optimizer_,NewOptimizer)
                obj.Optimizer_ = NewOptimizer;
            end
        end
        function val = get.Optimizer(obj)
            val = obj.Optimizer_;
        end
        function val = get.OptimizerParameters(obj)
            val = obj.OptimizerParameters_;
        end
        function obj = set.GradientThreshold(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','positive','nonnan'},'','GradientThreshold');
            obj.GradientThreshold = Value;
        end        
        function obj = set.GradientThresholdMethod(obj,Value)           
            validstr = validatestring(Value,{'l2norm','global-l2norm','absolute-value'},'','GradientThresholdMethod');
            obj.GradientThresholdMethod = string(validstr);
        end
        function obj = set.L2RegularizationFactor(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','nonnegative','finite'},'','L2RegularizationFactor');
            obj.L2RegularizationFactor = Value;
        end
        function obj = set.UseDevice(obj,Value)
            validatestring(Value,{'gpu','cpu'},'','UseDevice');
            obj.UseDevice = lower(string(Value));
        end
    end
end