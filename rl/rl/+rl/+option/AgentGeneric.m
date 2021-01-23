classdef AgentGeneric 
    % AGENTGENERIC: Creates agent options for all agent types
    
    % Copyright 2018 The MathWorks, Inc.

    properties
        % Sample time
        SampleTime
        
        % Discount factor to apply to future rewards during training
        DiscountFactor
    end
    
    properties (Access = protected)
        Parser        
    end    

    methods
        function obj = AgentGeneric(varargin)
            parser = rl.util.createInputParser;
            
            addParameter(parser,'SampleTime',1.0);
            % Parameters for reward discounting
            addParameter(parser,'DiscountFactor',0.99);

            parse(parser,varargin{:});
            obj.Parser = parser;
            obj.SampleTime = parser.Results.SampleTime;
            obj.DiscountFactor = parser.Results.DiscountFactor;
        end
        function obj = set.SampleTime(obj,Value)
            if isnumeric(Value) && isscalar(Value) && isreal(Value) && (Value>0 || Value==-1) && ~isinf(Value)
                obj.SampleTime = Value;                
            else
                error(message('rl:general:SampleTimeError'));
            end
        end
        function obj = set.DiscountFactor(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','positive','>',0,'<=',1},'','DiscountFactor');
            obj.DiscountFactor = Value;
        end
    end
    
end