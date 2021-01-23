classdef rlTD3AgentOptions < rl.option.AgentMemoryTarget
    % rlTD3AgentOptions: Creates agent options for TD3    
    
    % Copyright 2019 The MathWorks Inc.
    
    properties
        % Parameters for exploration model
        ExplorationModel
        % Parameters for target policy smoothing model
        TargetPolicySmoothModel
        % How many steps policy is updated
        PolicyUpdateFrequency
    end
    
    methods
        function obj = rlTD3AgentOptions(varargin)
            obj = obj@rl.option.AgentMemoryTarget(varargin{:});
            
            parser = obj.Parser; 
            addParameter(parser,'ExplorationModel',rl.option.GaussianActionNoise('Variance',0.1));
            addParameter(parser,'TargetPolicySmoothModel',rl.option.GaussianActionNoise('Variance',0.2,'LowerLimit',-0.5,'UpperLimit',0.5));
            addParameter(parser,'PolicyUpdateFrequency',2);

            parse(parser,varargin{:})            
            obj.Parser = parser;
            obj.ExplorationModel = parser.Results.ExplorationModel;
            obj.TargetPolicySmoothModel = parser.Results.TargetPolicySmoothModel;
            obj.PolicyUpdateFrequency = parser.Results.PolicyUpdateFrequency;
            obj.TargetUpdateFrequency = parser.Results.TargetUpdateFrequency;
            
            parser.KeepUnmatched = false;
            parse(parser,varargin{:}) 
        end
        function obj = set.ExplorationModel(obj,Value)
            validateattributes(Value,{'rl.option.OrnsteinUhlenbeckActionNoise','rl.option.GaussianActionNoise'},{'scalar'},'','ExplorationModel');
            obj.ExplorationModel = Value;
        end
        function obj = set.TargetPolicySmoothModel(obj,Value)
            validateattributes(Value,{'rl.option.GaussianActionNoise'},{'scalar'},'','TargetPolicySmoothModel');
            obj.TargetPolicySmoothModel = Value;
        end
        function obj = set.PolicyUpdateFrequency(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','nonnegative','integer'},'','PolicyUpdateFrequency');
            obj.PolicyUpdateFrequency = Value;
        end
    end
    
    methods (Access = protected)
        function DefaultTargetUpdateFrequency = getDefaultTargetUpdateFrequency(~)
            DefaultTargetUpdateFrequency = 2;
        end
        function DefaultTargetSmoothFactor = getDefaultTargetSmoothFactor(~)
            DefaultTargetSmoothFactor = 5e-3;
        end
    end
end