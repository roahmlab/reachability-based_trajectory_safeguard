classdef rlACAgentOptions < rl.option.AgentGeneric
    % rlACAgentOptions: Creates agent options for AC
    
    % Copyright 2017-2018 The MathWorks, Inc.
    
    properties
        % Number of steps to look-ahead
        NumStepsToLookAhead
    end
    properties
        % Entropy Weight
        EntropyLossWeight
    end
    
    methods
        function obj = rlACAgentOptions(varargin)
            obj = obj@rl.option.AgentGeneric(varargin{:});   
            
            parser = obj.Parser;  
            parser.KeepUnmatched = false;
            addParameter(parser,'NumStepsToLookAhead',1);
            addParameter(parser,'EntropyLossWeight',0);
            parse(parser,varargin{:});
            
            obj.Parser = parser;
            obj.NumStepsToLookAhead = parser.Results.NumStepsToLookAhead;
            obj.EntropyLossWeight = parser.Results.EntropyLossWeight;
        end
                
        function obj = set.NumStepsToLookAhead(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','integer','positive','finite'},'','NumStepsToLookAhead');
            obj.NumStepsToLookAhead = Value;
        end
        
        function obj = set.EntropyLossWeight(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','nonnegative','finite','<=',1},'','EntropyLossWeight');
            obj.EntropyLossWeight = Value;
        end
    end
end