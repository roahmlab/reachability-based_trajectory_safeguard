classdef BaseSimOptions
    % BASESIMOPTIONS

    % Copyright 2017-2018 The MathWorks Inc.    
    
    properties
        % Stop on simulation error if "on". Otherwise, attempt to continue
        % simulations if run in a loop.
        StopOnError
    end 
    properties (Access = protected)
        UseParallel_
        ParallelizationOptions_
    end
    properties (Dependent)        
        % parallel options
        UseParallel
        ParallelizationOptions
    end 
    
    methods
        function obj = BaseSimOptions(varargin)
            %% Parsing name-value pairs
            parser = inputParser;
            
            obj.addParseParams(parser);                      
            parse(parser,varargin{:});             
            
            %% Property assignments  

            obj.UseParallel = parser.Results.UseParallel;
            obj.ParallelizationOptions = parser.Results.ParallelizationOptions;
            obj.StopOnError = parser.Results.StopOnError;
        end
        function obj = set.UseParallel(obj,val)
            validateattributes(val,{'logical','numeric'},{'scalar','real'},'','UseParallel');
            obj.UseParallel_ = logical(val);
        end
        function val = get.UseParallel(obj)
            val = obj.UseParallel_;
        end
        function obj = set.ParallelizationOptions(obj,val)
            val = obj.validateParallelizationOptions(val);
            obj.ParallelizationOptions_ = val;
        end
        function val = get.ParallelizationOptions(obj)
            val = obj.ParallelizationOptions_;
        end
        function obj = set.StopOnError(obj,Value)
            obj.StopOnError = validatestring(Value,["on","off"],'','StopOnError');
        end 
    end
    methods
        function opts = getSimulationOptions(obj)
            opts = obj;
        end
    end
    methods (Access = protected,Static)
        function val = validateParallelizationOptions(val)
            validateattributes(val,{'rl.option.ParallelSimulation'},{'scalar'},'','ParallelizationOptions');
            % cast back to parallel sim
            val = rl.option.ParallelSimulation(val);
        end
        function opts = getDefaultParallelOptions()
            opts = rl.option.ParallelSimulation();
        end
    end
    methods (Access = protected)
        function addParseParams(obj,parser)
            % parallel parameters
            addParameter(parser,'UseParallel',false);
            addParameter(parser,'ParallelizationOptions',...
                obj.getDefaultParallelOptions());
            
            % stop on error parameter
            addParameter(parser,'StopOnError','on');
        end
    end
end