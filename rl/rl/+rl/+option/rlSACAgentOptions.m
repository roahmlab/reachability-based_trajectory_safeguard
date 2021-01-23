classdef rlSACAgentOptions < rl.option.AgentMemoryTarget
    % rlSACAgentOptions: Creates agent options for SAC    
    % Copyright 2019 The MathWorks Inc.
    
    properties
        % 
        EntropyWeightOptions
    end   
    methods
        function obj = rlSACAgentOptions(varargin)
            obj = obj@rl.option.AgentMemoryTarget(varargin{:});
            parser = obj.Parser; 
            addParameter(parser, 'EntropyWeightOptions',rl.option.EntropyWeightOptions);
            
            parse(parser,varargin{:});
            obj.Parser = parser;
            obj.EntropyWeightOptions = parser.Results.EntropyWeightOptions;
            
            parser.KeepUnmatched = false;
            parse(parser,varargin{:});
        end
        function obj = set.EntropyWeightOptions(obj,Value)
            validateattributes(Value,{'rl.option.EntropyWeightOptions'},{'scalar'},'','EntropyWeightOptions')
            obj.EntropyWeightOptions = Value;
        end
    end
end