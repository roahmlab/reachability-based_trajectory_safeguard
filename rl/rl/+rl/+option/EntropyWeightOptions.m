classdef EntropyWeightOptions
    % EntropyWeightOptions: Creates options for EntropyWeight    
    % Copyright 2019 The MathWorks Inc.
    
    properties
        % Initial EntropyWeight
        EntropyWeight
        
        % LearningRate
        LearnRate
        
        % desired target entropy of the policy
        TargetEntropy 
        
        % GradientThreshold
        GradientThreshold
    end
    
    properties (Dependent)
        Optimizer
        OptimizerParameters
    end

    properties (Access = private)
        Optimizer_
        OptimizerParameters_
    end
    
    methods
        function obj = EntropyWeightOptions(varargin)
            parser = inputParser;
            addParameter(parser,'EntropyWeight',1);
            addParameter(parser,'TargetEntropy',-1);
            addParameter(parser,'Optimizer',"adam");
            addParameter(parser,'LearnRate',3e-4);
            addParameter(parser,'GradientThreshold',Inf);
            
            parse(parser, varargin{:});
            
            obj.EntropyWeight = parser.Results.EntropyWeight;
            obj.LearnRate = parser.Results.LearnRate;
            obj.TargetEntropy = parser.Results.TargetEntropy;
            obj.GradientThreshold = parser.Results.GradientThreshold;
            obj.Optimizer = parser.Results.Optimizer;
        end
        function obj = set.EntropyWeight(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','positive','finite'},'','EntropyWeight');
            obj.EntropyWeight = Value;
        end
        function obj = set.LearnRate(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','nonnegative','finite'},'','LearnRate');
            obj.LearnRate = Value;
        end
        function obj = set.TargetEntropy(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','nonzero','finite'},'','TargetEntropy');
            if Value < 0 && Value ~= -1
                error(message('rl:agent:errEntropyNegative'));
            end
            obj.TargetEntropy = Value;
        end
        function obj = set.GradientThreshold(obj,Value)
             validateattributes(Value,{'numeric'},{'scalar','real','positive','nonnan'},'','GradientThreshold');
            obj.GradientThreshold = Value;
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
    end
end