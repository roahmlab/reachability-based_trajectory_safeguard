classdef ParallelTraining < rl.option.ParallelSimulation
% PARALLELTRAINING

% Revised: 10-26-2018
% Copyright 2018 The MathWorks, Inc.

    properties
        % Mode
        %   Parallelization mode:
        %       async
        %       sync
        Mode
        
        % DataToSendFromWorkers
        %   Type of data to send to the workers. Supported values include:
        %       Experiences
        %       Gradients
        DataToSendFromWorkers
        
        % StepsUntilDataIsSent
        %   Dictates when the worker sends data to the host and receives
        %   updated parameters. The size of the data send to the host will
        %   be [StepsUntilDataIsSent 1]. If -1, the worker will wait until
        %   the end of the episode to send data to the host.
        StepsUntilDataIsSent
    end
    properties (Hidden)
        % BlockWorkersUntilParametersAreReceived
        %   once workers send data back to host, block until parameters are
        %   received if true
        BlockWorkersUntilParametersAreReceived (1,1) logical = true
    end
    properties (Constant,Hidden)
        ValidParallelModes = {'sync','async'}
    end
    methods
        function this = ParallelTraining(varargin)
            this = this@rl.option.ParallelSimulation(varargin{:});
            
            parser = inputParser;
            this.addParseParams(parser);
            parse(parser,varargin{:});
            
            r = parser.Results;
            this.Mode = r.Mode;
            this.DataToSendFromWorkers = r.DataToSendFromWorkers;
            this.StepsUntilDataIsSent = r.StepsUntilDataIsSent;
        end
    end
    methods
        function this = set.Mode(this,val)
            val = validatestring(val,this.ValidParallelModes,'','Mode');
            this.Mode = string(val);
        end
        function this = set.DataToSendFromWorkers(this,val)
            validatestring(val,{'experiences','gradients'},'','DataToSendFromWorkers');
            this.DataToSendFromWorkers = string(val);
        end
        function this = set.StepsUntilDataIsSent(this,val)
            if ~(isscalar(val) && val == -1)
                validateattributes(val,{'numeric'},{'integer','scalar','positive'},'','StepsUntilDataIsSent');
            end
            this.StepsUntilDataIsSent = val;
        end
    end
    methods (Access = protected)
        function addParseParams(obj,parser)
            addParseParams@rl.option.ParallelSimulation(obj,parser)
            addParameter(parser,'Mode',"sync");
            addParameter(parser,'DataToSendFromWorkers',"experiences");
            addParameter(parser,'StepsUntilDataIsSent',-1);
        end
    end
end