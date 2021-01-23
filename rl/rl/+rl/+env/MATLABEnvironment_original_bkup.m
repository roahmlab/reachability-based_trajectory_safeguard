classdef MATLABEnvironment < rl.env.AbstractEnv
% MATLABENVIRONMENT

% Copyright 2018 The MathWorks, Inc.

    events
        % use this event to update other components (e.g. visualizers) due
        % to changes in the environment
        EnvUpdated
    end
    methods
        function validateEnvironment(this)
            % Validate environment by performing a short simulation.
            % The agent is reset, step one times with a random action from
            % ActionInfo spec, and reset again at the end.
            
            % validate output from reset function
            try
                InitialObservation = reset(this);
            catch ex
                me = MException(message('rl:env:errMatlabEnvResetCannotEvaluate'));
                throw(addCause(me,ex))
            end
            if ~isa(InitialObservation,'cell')
                InitialObservation = {InitialObservation};
            end
            if ~isDataValid(this.ObservationInfo,InitialObservation)
                error(message('rl:env:errMatlabEnvResetCheckObsInfoNotMatch'))
            end
            % validate output from step function, input a sample from
            % ActionInfo
            testAction = this.ActionInfo.usample;
            % unwrap cell for single action
            if numel(testAction) == 1
                testAction = testAction{:};
            end
            % throw error if step cannot be execute
            try
                Observation = step(this,testAction);
            catch ex
                me = MException(message('rl:env:errMatlabEnvStepCannotEvaluate'));
                throw(addCause(me,ex))
            end
            if ~isa(Observation,'cell')
                Observation = {Observation};
            end
            if ~isDataValid(this.ObservationInfo,Observation)
                error(message('rl:env:errMatlabEnvStepCheckObsInfoNotMatch'))
            end
            % clean up: reset the environment
            InitialObservation = reset(this); %#ok<NASGU>
        end
    end
    methods (Access = protected)
        function [experiences,varargout] = simWithPolicyImpl(env,policy,opts,varargin)
            % MATLAB simulation loop
            
            % initialize a cell array to store experiences per
            % episode
            N = opts.NumSimulations;
            expcell  = cell(N,1);
            siminfos = cell(N,1);
            
            usePCT = opts.UseParallel && matlab.internal.parallel.isPCTInstalled();
            % don't check worker id here g2036797
            if usePCT
                
                % distribute the "constant" data to the workers
                s.env = env;
                s.policy = policy;
                s.opts = opts;
                s.usePCT = opts.UseParallel && matlab.internal.parallel.isPCTInstalled();
                c = parallel.pool.Constant(s);
                
                simfh = @(simCount) simLoop(c.Value.env,c.Value.policy,c.Value.opts,simCount,c.Value.usePCT);
                
                % create the future objects with parfeval
                F = parallel.Future.empty(opts.NumSimulations,0);
                for simCount = 1:N
                    F(simCount) = parfeval(simfh,4,simCount); 
                end
                
                % tell the world that the sims have been distributed
                notifySimsDistributed(env);
                
                % fetch the outputs from parfeval
                for simCount = 1:N
                    if ~policy.TerminateSimulation
                        [~,expcell{simCount},epinfo,siminfos{simCount},workerID] = fetchNext(F);
                        % tell the world an episode has finished
                        notifyEpisodeFinished(env,epinfo,siminfos{simCount},simCount,workerID);
                    else
                        cancel(F)
                        break;
                    end
                end
            else
                % tell the world that the sims have been distributed
                notifySimsDistributed(env);
                for simCount = 1:N
                    % exit early
                    if policy.TerminateSimulation
                        break;
                    end
                    % run the sim loop
                    [expcell{simCount},epinfo,siminfos{simCount}] = simLoop(env,policy,opts,simCount,usePCT);
                    % tell the world an episode has finished
                    notifyEpisodeFinished(env,epinfo,siminfos{simCount},simCount);
                end
            end
            % put the experiences into an array
            experiences = expcell(1:simCount);
            
            if nargout > 1
                varargout{1} = vertcat(siminfos{:});
            end
        end
        function setupForSimImpl(env,opts)
            % check if options have changed and we need to re-setup
            if opts.UseParallel
                % grab the parallel pool
                ppool = gcp();
                
                % extract the parallel options
                paropts = opts.ParallelizationOptions;
                afiles = paropts.AttachedFiles;
                setupFcn = paropts.SetupFcn;
                xferbasews = paropts.TransferBaseWorkspaceVariables;
                
                % make sure this class is distributed to the workers as
                % well
                classfile = which(class(env));
                rclassfiles = matlab.codetools.requiredFilesAndProducts(classfile);
                
                if isempty(afiles)
                    afiles = {};
                else
                    afiles = cellstr(afiles);
                end
                files = [rclassfiles(:);afiles(:)];
                
                % Check for files to attach to the pool. Note files
                % attached to the pool are always using native paths,
                % convert any paths we have to native format
                files = regexprep(files,'[\\/]',filesep);
                nFiles = numel(files);
                if nFiles > 0
                    %Have files to attach
                    aFiles = ppool.AttachedFiles;
                    idx = ismember(files,aFiles);
                    if any(idx)
                        % Update files that already attached to the
                        % pool
                        updateAttachedFiles(ppool);
                    end
                    if any(~idx)
                        % Attach missing files to the pool
                        addAttachedFiles(ppool,files(~idx));
                    end
                end
                
                % get base ws variables
                if strcmp(xferbasews,"on")
                    baseWSVars = localGetBaseWSVars();
                else
                    baseWSVars = struct;
                end
                
                % Setup workers
                F = parfevalOnAll(@localSetupWorker,0,baseWSVars,setupFcn);
                
                % handle any errors
                try
                    fetchOutputs(F);
                catch ex
                    throwAsCaller(ex);
                end
            end
        end
        function cleanupForSimImpl(env)
            if env.SetupForSimOptions.UseParallel && ...
                    (rl.util.getWorkerID() == 0) && ...
                    matlab.internal.parallel.isPCTInstalled()
                % extract the parallel options
                cleanupFcn = env.SetupForSimOptions.ParallelizationOptions.CleanupFcn;
                
                % cleanup workers
                F = parfevalOnAll(@localCleanupWorker,0,cleanupFcn);
                
                % handle any errors
                try
                    fetchOutputs(F);
                catch ex
                    throwAsCaller(ex);
                end
            end
        end
        function notifyEnvUpdated(this)
            % call this function any time you want to tell the world that
            % the environment has been updated
            
            % fire EnvUpdated
            notify(this,'EnvUpdated');
            
            % fire user-defined callback for environments that do plotting
            % internally
            envUpdatedCallback(this);
        end
        function envUpdatedCallback(this) %#ok<MANU>
            % overload this function to execute code once the environment
            % has been updated
        end
    end
    methods (Access = private)
        function [experiences,episodeInfo,simInfo,workerID] = simLoop(env,policy,opts,simCount,usePCT)
            % runs the simulation loop
                        
            % do stuff on the policy pre sim
            preSimFcn(policy,simCount);
            
            simInfo = struct('SimulationError',[]);
            try
                % reset the environment
                observation = reset(env);

                % get the initial action from the initial observation
                action = getInitialAction(policy,observation);

                % run the simulations
                for stepCount = 1:opts.MaxSteps

                    % exit the simulation if termination requested
                    if policy.TerminateSimulation
                        break;
                    end

                    % step the environment
                    [observation,reward,isdone,loggedSignals] = step(env,action); %#ok<ASGLU>

                    % step the policy
                    action = step(policy,observation,reward,isdone);

                    if isdone
                        break;
                    end
                end
            catch ex
                if strcmpi(opts.StopOnError,"off")
                    simInfo.SimulationError = ex;
                else
                    rethrow(ex);
                end
            end
            
            % do stuff on the policy post sim
            postSimFcn(policy,simCount,simInfo);
            
            % collect the experiences (will return empty if logger was
            % not attached)
            experiences = getExperiences(policy);
            
            % collect the episode info
            episodeInfo = getEpisodeInfo(policy);
            
            % get the workerID
            workerID = rl.util.getWorkerID(usePCT);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Local Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localSetupWorker(vars,setupFcn)
% setup workers for parallel training

% REVISIT a lot of this was a hold over from simulink setup for frestimate.
% We should only have to run the setupFcn

% first clean up the workers
evalin('base','clear');

% distribute base ws vars
addedvars = localDistributeToBaseWS(vars);

% build the simWorker struct
workerData.varsAdded = addedvars;

% install the sim worker
localSetGetWorkerData(workerData);

% run the user setup fcn
if ~isempty(setupFcn)
    setupFcn();
end
end

function localCleanupWorker(cleanupFcn)
% cleanup the workers after training

% run the user cleanup fcn
if ~isempty(cleanupFcn)
    cleanupFcn();
end

% get the sim worker data that was created in parSetupWorker
workerData = localSetGetWorkerData();

% Restore any changes made during setupWorkers

% Clear variables added to workspace
if ~isempty(workerData.varsAdded)
    cmd = sprintf('%s ',workerData.varsAdded{:});
    cmd = sprintf('clear %s', cmd);
    evalin('base',cmd)
end
end
function s = localSetGetWorkerData(s)
% install simWorker globally
persistent data_;
if nargin
    data_ = s;
end
s = data_;
end
function vars = localGetBaseWSVars()
% get the base ws variables
varList = evalin('base','who');

vars = struct;
for i = 1 : numel(varList)
    val = evalin('base',varList{i});
    if ~localIsNonSerializable(val)
        vars.(varList{i}) = val;
    end
end
end
function addedvars = localDistributeToBaseWS(vars)
% check for vars already in the base ws (may be opened by the model)
wsvars = evalin('base','who');
fs = fields(vars);
addedvars = setdiff(fs,wsvars);

for i = 1:numel(addedvars)
    varname = addedvars{i};
    varval  = vars.(varname);
    assignin('base',varname,varval);
end
end
function val = localIsNonSerializable(var)
% deternine if a variable is serializeable
val = ...
    isa(var,'Composite') || ...
    isa(var,'parallel.internal.customattr.CustomPropTypes') || ...
    isa(var,'distributed') || ...
    isa(var,'qeWorkingDir');
end