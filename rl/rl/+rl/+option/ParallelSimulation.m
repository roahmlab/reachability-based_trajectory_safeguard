classdef ParallelSimulation
% PARALLELSIMULATION

% Revised: 3-20-2019
% Copyright 2018 The MathWorks, Inc.

    properties
        % REVISIT: in the async case, the random seeds may not help as the
        % execution on each worker is not deterministic, e.g. worker 1 may
        % finish simulation 1 before worker 2 in scenario 1, however this
        % may not be the case in scenario 2, resulting in the host
        % receiving experiences/gradients in a different order and the host
        % sending back parameters that aren't the same as in scenario 1
        
        % WorkerRandomSeeds
        %   -2              : Do not put a random seed on the worker
        %   -1              : Put a unique random seed on each worker 
        %                     corresponding to worker id
        %   [1 2 3 ... N]   : Manually specify the random seed for each
        %                     worker. The number of elements MUST match the
        %                     number of workers
        WorkerRandomSeeds
        
        % TransferBaseWorkspaceVariables
        %   When 'on' variables used in models and defined in the base
        %   workspace will be transferred to the parallel workers
        TransferBaseWorkspaceVariables
        
        % AttachedFiles
        %   Additional files to attach to the parallel pool
        AttachedFiles
        
        % SetupFcn
        %   A function handle to run once per worker before training
        %   starts. The function handle must not require any input
        %   arguments.
        SetupFcn
        
        % CleanupFcn
        %   A function handle to run once per worker after training
        %   ends. The function handle must not require any input
        %   arguments.
        CleanupFcn
    end
    methods
        function this = ParallelSimulation(varargin)
            
            if nargin == 1 && isa(varargin{1},'rl.option.ParallelSimulation')
                that = varargin{1};
                this = rl.option.ParallelSimulation(...
                    'WorkerRandomSeeds',that.WorkerRandomSeeds,...
                    'AttachedFiles',that.AttachedFiles,...
                    'SetupFcn',that.SetupFcn,...
                    'CleanupFcn',that.CleanupFcn,...
                    'TransferBaseWorkspaceVariables',that.TransferBaseWorkspaceVariables);
            else
            
                parser = inputParser;
                this.addParseParams(parser);
                parse(parser,varargin{:});

                r = parser.Results;
                this.WorkerRandomSeeds = r.WorkerRandomSeeds;

                this.AttachedFiles = r.AttachedFiles;
                this.SetupFcn = r.SetupFcn;
                this.CleanupFcn = r.CleanupFcn;
                this.TransferBaseWorkspaceVariables = r.TransferBaseWorkspaceVariables;
            end
        end
    end
    methods
        function this = set.WorkerRandomSeeds(this,val)
            if ~(numel(val) == 1 && (val == -1 || val == -2))
                validateattributes(val,{'numeric'},{'integer','vector','nonnegative'},'','WorkerRandomSeeds');
            end
            this.WorkerRandomSeeds = val;
        end
        function this = set.AttachedFiles(this,val)
            if ~isempty(val)
                if iscellstr(val) %#ok<ISCLSTR>
                    val = string(val);
                end
                validateattributes(val,{'string'},{},'','AttachedFiles');
            end
            this.AttachedFiles = val;
        end
        function this = set.SetupFcn(this,val)
            if ~isempty(val)
                validateattributes(val,{'function_handle'},{'scalar'},'','SetupFcn');
                if nargin(val)
                    error(message('rl:general:ParallelTrainSetupFcnInvalidIputs'));
                end
            end
            this.SetupFcn = val;
        end
        function this = set.CleanupFcn(this,val)
            if ~isempty(val)
                validateattributes(val,{'function_handle'},{'scalar'},'','CleanupFcn');
                if nargin(val)
                    error(message('rl:general:ParallelTrainCleanupFcnInvalidIputs'));
                end
            end
            this.CleanupFcn = val;
        end
        function this = set.TransferBaseWorkspaceVariables(this,val)
            validatestring(val,["on","off"],'','TransferBaseWorkspaceVariables');
            this.TransferBaseWorkspaceVariables = string(val);
        end
    end
    methods (Access = protected)
        function addParseParams(obj,parser) %#ok<INUSL>
            addParameter(parser,'WorkerRandomSeeds',-1);
            addParameter(parser,'AttachedFiles',[]);
            addParameter(parser,'SetupFcn',[]);
            addParameter(parser,'CleanupFcn',[]);
            addParameter(parser,'TransferBaseWorkspaceVariables',"on");
        end
    end
end