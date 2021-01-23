classdef rlDDPGAgentOptions < rl.option.AgentMemoryTarget
    % rlDDPGAgentOptions: Creates agent options for DDPG    
    
    % Copyright 2017-2018 The MathWorks Inc.
    
    properties
        % Parameters for Ornstein Uhlenbeck Action Noise
        NoiseOptions
    end
    methods
        function obj = rlDDPGAgentOptions(varargin)
            obj = obj@rl.option.AgentMemoryTarget(varargin{:});
            
            parser = obj.Parser; 
            addParameter(parser,'NoiseOptions',rl.option.OrnsteinUhlenbeckActionNoise);

            parse(parser,varargin{:})            
            obj.Parser = parser;
            obj.NoiseOptions = parser.Results.NoiseOptions; 
            
            parser.KeepUnmatched = false;
            parse(parser,varargin{:}) 
        end
        function obj = set.NoiseOptions(obj,Value)
            validateattributes(Value,{'rl.option.OrnsteinUhlenbeckActionNoise','rl.option.GaussianActionNoise'},{'scalar'},'','NoiseOptions');
            obj.NoiseOptions = Value;
        end             
    end
end