classdef TaskSpec < handle & matlab.mixin.Heterogeneous
% TASKSPEC
%
% This class is responsible for building tasks and running them either on
% the host machine or on a thread via some parallelization mechanism

% Revised: 6-25-2019
% Copyright 2019 The MathWorks, Inc.

    properties (Transient)
        % Allow the task to run on a worker. Set to false if you want to
        % force the task to run on the host.
        RunOnWorker (1,1) logical = true
    end
    properties (SetAccess = private,Transient)
        % outputs from running the task
        Outputs cell = {};
        % is the task finished
        FinishedTask (1,1) logical = false
        % parfeval futures for worker task runs
        Future = []
    end
    properties (Access = private,Transient)
        % parallel pool object
        Pool = []
    end
    methods 
        function delete(this)
            % make sure the futures are cleaned up
            delete(this.Future);
            
            % restore the state of the feval callback queue
            localResetFevalQueueCallbacks(this.Pool);
        end
    end
    methods (Sealed)
        function run(this)
            % run an array of task specs
            
            % put the tasks in a row vector for convenience
            this = this(:)';
            % get the required num workers for each task
            nw = arrayfun(@getNumRequiredWorkers,this);
            
            wtaskIdx = nw > 0 & [this.RunOnWorker];
            htaskIdx = ~wtaskIdx;
            
            % attach the pool to the tasks that will be run on workers
            if matlab.internal.parallel.isPCTInstalled()
                p = gcp('nocreate');
            else
                p = [];
            end
            [this(wtaskIdx).Pool] = deal(p);
            
            % flush the state of the feval callback queue. We manage the
            % feval callback queue to make sure afterEach functions run as
            % expected since the parsim manager will turn the callback
            % queue off if the simulations are not run in background
            localEnableFevalQueueCallbacks(p,[]);
            
            % reorder the tasks:
            %   1. specified # of required workers (e.g. 1, 2 ...)
            %   2. host (e.g. 0)
            this = [this(wtaskIdx),this(htaskIdx)];
            
            % run the tasks in order
            for task = this
                runScalarTask(task);
            end
        end
        function spmdrun(this)
            % run an array of task specs using spmd rather than parfeval. 
            
            % put the tasks in a row vector for convenience
            this = this(:)';
            
            % get the required num workers for each task
            nw = arrayfun(@getNumRequiredWorkers,this);
            
            wtaskIdx = nw > 0 & [this.RunOnWorker];
            
            assert(all(wtaskIdx),...
                'spmdrun can only be used on TaskSpecs that are configured to run in parallel');
            
            N = numel(this);
            
            % create composites for each task
            composite_task = Composite();
            for i = 1:N
                composite_task{i} = this(i);
            end
            
            % start the SPMD block
            spmd (N)
                % turn off dead lock detection REVISIT can we query the
                % initial state of this feature to reset at the end of the
                % spmd block?
                % mpiSettings('DeadlockDetection','off');
                
                
                runDirect(composite_task);
                
                % flush any remaining messages to avoid spmd warnings
                labBarrier
                while labProbe
                    data2discard = labReceive(); %#ok<NASGU>
                end
                
                % capture the outputs from each task
                outs = composite_task.Outputs;
            end
            
            % spmd is blocking so mark the tasks as finished once we get
            % here
            [this.FinishedTask] = deal(true);
            
            % attach the outputs
            [this.Outputs] = outs{:};
        end
        function waitForAllTasks2Finish(this)
            % wait until all the tasks report they are finished
            futures = [this.Future];
            while ~all([this.FinishedTask])
                if ~isempty(futures)
                    err = [futures.Error];
                    idx = ~isempty(err);
                    if any(idx)
                        e = err(idx);
                        throw(e(1));
                    end
                end
                pause(0.1);
            end
        end
    end
    methods (Access = protected,Abstract)
        % build the task
        task = buildTask(this)
        % get the required number of workers for the task:
        %   0: no workers required (run on host)
        %   1: 1 worker required
        % Based on this number, the tasks will be ordered appropriately
        n = getNumRequiredWorkers(this)
        % get the number of outputs for the task
        n = getNumOutputs(this)
    end
    methods (Access = ?qe.rl.TaskTestCase)
        % helper testing methods
        function n = hGetNumOutputs(this)
            n = getNumOutputs(this);
        end
    end
    methods (Access = protected)
        function varargout = internal_run(this)
            % build and runs a scalar task
            task = buildTask(this);
            clnup = onCleanup(@() delete(task));
            [varargout{1:nargout}] = run(task);
        end
        function runDirect(this)
            % build and run the task on the current executing thread
            [this.Outputs{1:getNumOutputs(this)}] = internal_run(this);
            this.FinishedTask = true;
        end
        function runOnWorker(this,varargin)
            % build and run the task on a worker
                       
            % make sure the feval callback queue is turned on
            localEnableFevalQueueCallbacks(this.Pool,true);
            
            this.Future = parfeval(this.Pool,@internal_run,getNumOutputs(this),this);
            afterEach(this.Future,@(varargin) afterEachFutureAssign(this,varargin{:}),0);
        end
        function afterEachFutureAssign(this,varargin)
            % fetch the output of the future and assign it to Outputs. Mark
            % the task as done
            this.Outputs = varargin;
            this.FinishedTask = true;
        end
        function runScalarTask(this)
            % run a single task
            this.FinishedTask = false;
            % flush the outputs
            this.Outputs = {};
            if isempty(this.Pool) || (getNumRequiredWorkers(this) == 0)
                runDirect(this);
            else
                runOnWorker(this);
            end
        end
    end
end

%% local functions for managing the parfeval callback queue
function oldval = localEnableFevalQueueCallbacks(p,val)
% enable/disable feval queue callbacks. Needed for afterEach on the
% parfeval futures since the Simulink engine will shut disable the callback
% queue if not run in background mode
persistent oldval_
if nargin 
    if isempty(val)
        % enable resetting of oldval to []
        oldval_ = val;
    elseif ~isempty(p) && (isempty(oldval_) || val ~= oldval_)
        % turn on/off the feval queue callbacks, keep the oldval in the
        % persistent variable
        oldval_ = hToggleCallbacks(p.FevalQueue,val);
    end
end
oldval = oldval_;
end
function localResetFevalQueueCallbacks(p)
% reset the feval queue callback state to its original value
oldval = localEnableFevalQueueCallbacks();
localEnableFevalQueueCallbacks(p,oldval);
end