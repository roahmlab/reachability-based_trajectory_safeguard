classdef AgentMemoryTarget < rl.option.AgentGeneric
    % AGENTMEMORYTARGET: Creates agent options for agents with memory and targets
    
    % Copyright 2018 The MathWorks, Inc.

    properties
        % Smoothing factor that determines how a target model used for
        % prediction is updated with weights from the prediction model used
        % for training. Set it to 1 if the weights of the prediction model
        % need to be copied over to the target model
        TargetSmoothFactor
        
        % Number of steps after which the target model is updated
        TargetUpdateFrequency
                
        % Reset experience buffer before training
        ResetExperienceBufferBeforeTraining
        
        % Save experience buffer with the agent when calling save
        SaveExperienceBufferWithAgent = false
    end
    
    properties (Dependent,Hidden)
        % Strategy to use to update the target network during training
        TargetUpdateMethod % Kept for backward compatibility
    end
    
    properties (Dependent)
        % Batch size for agent training
        MiniBatchSize
        
        % Number of steps to look-ahead
        NumStepsToLookAhead
        
        % Experience buffer length
        ExperienceBufferLength
    end
    
    properties (Access = protected)
        % Batch size for agent training
        MiniBatchSize_
        
        % Number of steps to look-ahead
        NumStepsToLookAhead_
        
        % Experience buffer length
        ExperienceBufferLength_
    end
    methods
        function obj = AgentMemoryTarget(varargin)
            obj = obj@rl.option.AgentGeneric(varargin{:});
            
            parser = obj.Parser;            
            addParameter(parser,'MiniBatchSize',64);            
            addParameter(parser,'TargetSmoothFactor',obj.getDefaultTargetSmoothFactor);
            addParameter(parser,'TargetUpdateFrequency',obj.getDefaultTargetUpdateFrequency);
            addParameter(parser,'ResetExperienceBufferBeforeTraining',true);
            addParameter(parser,'SaveExperienceBufferWithAgent',false);
            addParameter(parser,'NumStepsToLookAhead',1);
            addParameter(parser,'ExperienceBufferLength',10000);
            addParameter(parser,'TargetUpdateMethod',[]);
            parse(parser,varargin{:});
            
            obj.Parser = parser;
            obj.NumStepsToLookAhead = parser.Results.NumStepsToLookAhead;  
            obj.MiniBatchSize = parser.Results.MiniBatchSize;            
            obj.TargetSmoothFactor = parser.Results.TargetSmoothFactor;
            obj.TargetUpdateFrequency = parser.Results.TargetUpdateFrequency;
            obj.ResetExperienceBufferBeforeTraining = parser.Results.ResetExperienceBufferBeforeTraining;
            obj.SaveExperienceBufferWithAgent = parser.Results.SaveExperienceBufferWithAgent;
            obj.ExperienceBufferLength = parser.Results.ExperienceBufferLength;
            if ~isempty(parser.Results.TargetUpdateMethod)
                obj.TargetUpdateMethod = parser.Results.TargetUpdateMethod;
            end
        end
        function obj = set.TargetSmoothFactor(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','positive','<=',1},'','TargetSmoothFactor');
            obj.TargetSmoothFactor = Value;
        end
        function obj = set.TargetUpdateFrequency(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','nonnegative','integer'},'','TargetUpdateFrequency');
            obj.TargetUpdateFrequency = Value;
        end
        function obj = set.TargetUpdateMethod(obj,Value)
            Criteria = rl.option.AgentMemoryTarget.getTargetUpdateMethod;
            validatestring(Value,Criteria,'','TargetUpdateMethod');
            switch Value
                case "smoothing"
                    % force update at every step
                    obj.TargetUpdateFrequency = 1;
                case "periodic"
                    % prevent soft update
                    obj.TargetSmoothFactor = 1;
                case "periodic-smoothing"
                    % nothing to do here
                case "copy"
                    % force update at every step and prevent soft update
                    obj.TargetUpdateFrequency = 1;
                    obj.TargetSmoothFactor = 1;
            end
        end
        function obj = set.ResetExperienceBufferBeforeTraining(obj,Value)
            validateattributes(Value,{'logical','numeric'},{'scalar','real','finite'},'','ResetExperienceBufferBeforeTraining');
            obj.ResetExperienceBufferBeforeTraining = Value;
        end
        
        function obj = set.SaveExperienceBufferWithAgent(obj,Value)
            validateattributes(Value,{'logical','numeric'},{'scalar','real','finite'},'','SaveExperienceBufferWithAgent');
            obj.SaveExperienceBufferWithAgent = Value;
        end
        
        function obj = set.MiniBatchSize(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','integer','positive','finite'},'','MiniBatchSize');
            if Value > obj.ExperienceBufferLength_
                error(message('rl:general:errBufferLengthLessThanMiniBatchSize'));
            end
            obj.MiniBatchSize_ = Value;
        end
        
        function obj = set.NumStepsToLookAhead(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','integer','positive','finite'},'','NumStepsToLookAhead');
            if Value > obj.ExperienceBufferLength_
                error(message('rl:general:AgentOptionsErrorBufferLengthLessThanNStep'));
            end
            obj.NumStepsToLookAhead_ = Value;
        end
        
        function obj = set.ExperienceBufferLength(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','integer','positive'},'','ExperienceBufferLength');
            if Value < obj.NumStepsToLookAhead_
                error(message('rl:general:AgentOptionsErrorBufferLengthLessThanNStep'));
            end
            if Value < obj.MiniBatchSize_
                error(message('rl:general:errBufferLengthLessThanMiniBatchSize'));
            end
            obj.ExperienceBufferLength_ = Value;
        end
        
        function Options = get.MiniBatchSize(this)
            Options = this.MiniBatchSize_;
        end
        
        function Options = get.NumStepsToLookAhead(this)
            Options = this.NumStepsToLookAhead_;
        end
        
        function Options = get.ExperienceBufferLength(this)
            Options = this.ExperienceBufferLength_;
        end
        function val = get.TargetUpdateMethod(this)
            freq = this.TargetUpdateFrequency;
            fact = this.TargetSmoothFactor;
            if (freq > 1) && (fact < 1)
                val = "periodic-smoothing";
            elseif (freq == 1) && (fact < 1)
                val = "smoothing";
            elseif (freq > 1) && (fact == 1)
                val = "periodic";
            else
                val = "copy";
            end
        end
    end
    
    methods (Static, Access = private)
        function Method = getTargetUpdateMethod()
            Method = ["smoothing","periodic","periodic-smoothing","copy"];
        end
    end
    methods (Access = protected)
        % NOTE: overwrite if need to change the default behavior of subclass.
        function DefaultTargetUpdateFrequency = getDefaultTargetUpdateFrequency(~)
            % E.g. TD3 changes to 2;
            DefaultTargetUpdateFrequency = 1;
        end
        function DefaultTargetSmoothFactor = getDefaultTargetSmoothFactor(~)
            % E.g. TD3 changes to 5e-3;
            DefaultTargetSmoothFactor = 1e-3;
        end
    end
end
