classdef rlPGAgentOptions < rl.option.AgentGeneric
    % rlPGAgentOptions: Creates agent options for PG
    
    % Copyright 2017-2018 The MathWorks Inc.
    
    properties
        % Whether uses baseline for learning
        UseBaseline
    end
    properties
        % Entropy Weight
        EntropyLossWeight
    end
    
    methods
        function obj = rlPGAgentOptions(varargin)
            obj = obj@rl.option.AgentGeneric(varargin{:});
            
            parser = obj.Parser;
            addParameter(parser,'UseBaseline',true);
            addParameter(parser,'EntropyLossWeight',0);
            parse(parser,varargin{:});                                  
            
            obj.Parser = parser;
            obj.UseBaseline = parser.Results.UseBaseline;
            obj.EntropyLossWeight = parser.Results.EntropyLossWeight;
            
            parser.KeepUnmatched = false;
            parse(parser,varargin{:})               
        end
        
        function obj = set.UseBaseline(obj,Value)
            validateattributes(Value,{'logical','numeric'},{'scalar','real'},'','UseBaseline');
            obj.UseBaseline = logical(Value);
        end
        
        function obj = set.EntropyLossWeight(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','real','nonnegative','finite','<=',1},'','EntropyLossWeight');
            obj.EntropyLossWeight = Value;
        end
    end
end