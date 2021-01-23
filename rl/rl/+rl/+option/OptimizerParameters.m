classdef OptimizerParameters
    % OptimizerParameters: Creates parameters for optimizer.
    
    % Copyright 2017-2018 The MathWorks Inc.
    
    properties (Dependent)
        Momentum
        Epsilon
        GradientDecayFactor
        SquaredGradientDecayFactor
    end
    properties (Access=private)
        Optimizer = 'adam';   
        
        Momentum_
        Epsilon_
        GradientDecayFactor_
        SquaredGradientDecayFactor_
    end
    methods
        function obj = OptimizerParameters(Optimizer)
            obj = setOptimizer(obj,Optimizer);
        end
        function obj = set.Momentum(obj,Value)
            if strcmpi(obj.Optimizer,{'sgdm'})
                validateattributes(Value,{'numeric'},{'real','nonnegative','<=',1},'','Momentum');
                obj.Momentum_ = Value;
            else
                obj.Momentum_ = "Not applicable";
            end
        end
        function obj = set.Epsilon(obj,Value)
            if any(strcmpi(obj.Optimizer,{'rmsprop','adam'}))
                validateattributes(Value,{'numeric'},{'scalar','real','positive'},'','Epsilon');
                obj.Epsilon_ = Value;
            else
                obj.Epsilon_ = "Not applicable";               
            end            
        end
        function obj = set.GradientDecayFactor(obj,Value)
            if any(strcmpi(obj.Optimizer,{'adam'}))
                validateattributes(Value,{'numeric'},{'real','scalar','nonnegative','<=',1},'','GradientDecayFactor');
                obj.GradientDecayFactor_ = Value;
            else
                obj.GradientDecayFactor_ = "Not applicable";
            end
        end        
        function obj = set.SquaredGradientDecayFactor(obj,Value)
            if any(strcmpi(obj.Optimizer,{'rmsprop','adam'}))
                validateattributes(Value,{'numeric'},{'real','nonnegative','<=',1},'','SquaredGradientDecayFactor');
                obj.SquaredGradientDecayFactor_ = Value;
            else
                obj.SquaredGradientDecayFactor_ = "Not applicable";
            end
        end 
        function val = get.Momentum(obj)
            val = obj.Momentum_;
        end
        function val = get.Epsilon(obj)
            val = obj.Epsilon_;         
        end
        function val = get.GradientDecayFactor(obj)
            val = obj.GradientDecayFactor_;
        end        
        function val = get.SquaredGradientDecayFactor(obj)
            val = obj.SquaredGradientDecayFactor_;
        end
        function obj = setOptimizer(obj,Optimizer)
            validstr = validatestring(Optimizer,{'adam','sgdm','rmsprop'},'','Optimizer');
            obj.Optimizer = validstr;
            switch validstr
                case {'adam'}
                    if strcmpi(obj.GradientDecayFactor,"Not applicable") || isempty(obj.GradientDecayFactor)
                        obj.GradientDecayFactor = 0.9;
                    end
                    if strcmpi(obj.SquaredGradientDecayFactor,"Not applicable") || isempty(obj.SquaredGradientDecayFactor)
                        obj.SquaredGradientDecayFactor = 0.999;
                    end                    
                    if strcmpi(obj.Epsilon,"Not applicable") || isempty(obj.Epsilon)
                        obj.Epsilon = 1e-8;
                    end                    
                    obj.Momentum = "Not applicable";
                case 'rmsprop'
                    if strcmpi(obj.SquaredGradientDecayFactor,"Not applicable") || isempty(obj.SquaredGradientDecayFactor)
                        obj.SquaredGradientDecayFactor = 0.999;
                    end                    
                    if strcmpi(obj.Epsilon,"Not applicable") || isempty(obj.Epsilon)
                        obj.Epsilon = 1e-8;
                    end
                    obj.GradientDecayFactor = "Not applicable";
                    obj.Momentum = "Not applicable";
                case 'sgdm'
                    if strcmpi(obj.Momentum,"Not applicable") || isempty(obj.Momentum)
                        obj.Momentum = 0.9;
                    end
                    obj.Epsilon = "Not applicable";
                    obj.GradientDecayFactor = "Not applicable";
                    obj.SquaredGradientDecayFactor = "Not applicable";
            end                                           
        end
        function optimizer = getOptimizer(obj)
            optimizer = obj.Optimizer;            
        end
    end
end