classdef rlSimulationOptions < rl.option.BaseSimOptions
    % RLSIMULATIONOPTIONS: Creates simulation options for RL agent.

    % Copyright 2017-2018 The MathWorks Inc.    
    
    properties
        % Number of steps for simulation
        MaxSteps
        
        % Number of simulations to run
        NumSimulations
    end
    methods
        function obj = rlSimulationOptions(varargin)
            
            obj = obj@rl.option.BaseSimOptions(varargin{:});
            
            parser = inputParser;
            obj.addParseParams(parser);
            parse(parser,varargin{:});                         
            
            obj.MaxSteps = parser.Results.MaxSteps;
            obj.NumSimulations = parser.Results.NumSimulations;
        end
        function obj = set.MaxSteps(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','integer','positive'},'','MaxSteps');
            obj.MaxSteps = Value;
        end 
        function obj = set.NumSimulations(obj,Value)
            validateattributes(Value,{'numeric'},{'scalar','integer','positive'},'','NumSimulations');
            obj.NumSimulations = Value;
        end 
    end
    methods (Access = protected)
        function addParseParams(obj,parser) 
            addParseParams@rl.option.BaseSimOptions(obj,parser);
            addParameter(parser,'MaxSteps',500);
            addParameter(parser,'NumSimulations',1);
        end
    end
end