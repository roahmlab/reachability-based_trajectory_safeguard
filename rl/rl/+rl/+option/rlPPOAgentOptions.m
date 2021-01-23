classdef rlPPOAgentOptions < rl.option.AgentGeneric
    % rlPPOAgentOptions: Creates agent options for PPO
    
    % Copyright 2019 The MathWorks, Inc.
    
    properties (Dependent)
        % Number of steps the agent interacts with the environment before 
        % learns from experience.
        ExperienceHorizon
        % The size of the mini-batch used for each training iteration.
        MiniBatchSize
    end
    properties (Access = private)
        % Number of steps the agent interacts with the environment before 
        % learns from experience.
        ExperienceHorizon_
        % The size of the mini-batch used for each training iteration.
        MiniBatchSize_
    end
    properties (Hidden)
        % Objective function for PPO agent to optimize. 
        % Currently supports only 'clipped' objective.
        % Place holder for fixed-kl, adaptive-kl objective.
        PPOObjective
        % 'clipped': limit the change in each policy update step by clipping
        % the action probability ratio between the new and old policies.
    end
    properties
        % Clip factor to limit the change in each policy update step.
        ClipFactor
        % Weight for entropy loss to promote policy exploration.
        EntropyLossWeight
        % Number of epochs. An epoch is the full pass of the training 
        % algorithm over the entire training set.
        NumEpoch
        % Method used for estimating advantage values for learning.
        AdvantageEstimateMethod
        % Smooth factor of the generalized advantage estimator.
        GAEFactor
    end
    
    
    methods
        function obj = rlPPOAgentOptions(varargin)
            obj = obj@rl.option.AgentGeneric(varargin{:});   
            
            parser = obj.Parser;  
            parser.KeepUnmatched = false;
            addParameter(parser,'ExperienceHorizon',512);
            addParameter(parser,'ClipFactor',0.2);
            addParameter(parser,'EntropyLossWeight',0.01);
            addParameter(parser,'AdvantageEstimateMethod','gae');
            addParameter(parser,'GAEFactor',0.95);
            addParameter(parser,'MiniBatchSize',128);
            addParameter(parser,'NumEpoch',3);
            parse(parser,varargin{:});
            
            obj.Parser = parser;
            obj.ExperienceHorizon = parser.Results.ExperienceHorizon;
            obj.ClipFactor = parser.Results.ClipFactor;
            obj.EntropyLossWeight = parser.Results.EntropyLossWeight;
            obj.AdvantageEstimateMethod = parser.Results.AdvantageEstimateMethod;
            obj.GAEFactor = parser.Results.GAEFactor;
            % TODO: support fixed-kl, adaptive-kl objectives
            obj.PPOObjective = 'clipped';
            obj.MiniBatchSize = parser.Results.MiniBatchSize;
            obj.NumEpoch = parser.Results.NumEpoch;
        end
        
        function obj = set.PPOObjective(obj,Value)
            validatestring(Value,{'clipped'},'','PPOObjective');
            obj.PPOObjective = lower(string(Value));
        end
        
        function obj = set.ExperienceHorizon(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','integer','positive','finite'},'','ExperienceHorizon');
            if Value < obj.MiniBatchSize_
                error(message('rl:agent:errPPOMiniBatchSizeGTExpHorizon'));
            end
            obj.ExperienceHorizon_ = Value;
        end
        
        function obj = set.ClipFactor(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','positive','finite','<',1},'','ClipFactor');
            obj.ClipFactor = Value;
        end
        
        function obj = set.EntropyLossWeight(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','nonnegative','finite','<=',1},'','EntropyLossWeight');
            obj.EntropyLossWeight = Value;
        end
        
        function obj = set.AdvantageEstimateMethod(obj,Value)
            validatestring(Value,{'gae','finite-horizon'},'','AdvantageEstimateMethod');
            obj.AdvantageEstimateMethod = lower(string(Value));
        end
        
        function obj = set.GAEFactor(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','nonnegative','finite','<=',1},'','GAEFactor');
            obj.GAEFactor = Value;
        end
        
        function obj = set.MiniBatchSize(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','integer','positive','finite'},'','MiniBatchSize');
            if Value > obj.ExperienceHorizon_
                error(message('rl:agent:errPPOMiniBatchSizeGTExpHorizon'));
            end
            obj.MiniBatchSize_ = Value;
        end
        
        function obj = set.NumEpoch(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','integer','positive'},'','NumEpoch');
            obj.NumEpoch = Value;
        end
        
        function Options = get.MiniBatchSize(this)
            Options = this.MiniBatchSize_;
        end
        
        function Options = get.ExperienceHorizon(this)
            Options = this.ExperienceHorizon_;
        end
    end
end