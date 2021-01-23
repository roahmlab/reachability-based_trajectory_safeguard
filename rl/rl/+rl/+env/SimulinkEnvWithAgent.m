classdef SimulinkEnvWithAgent < rl.env.AbstractEnv
    % SIMULINKENVWITHAGENT: This class defines a Simulink environment
    % with an agent block in the loop, representing the integrated RL
    % workflow
    
    % Copyright 2017-2018 The MathWorks Inc.
    
    properties (SetAccess = private)
        % ENV.MODEL
        %   Read only model name.
        %
        %   mdl = env.Model;
        Model string
        
        % ENV.AGENTBLOCK
        %   Read only path to the agent block.
        %
        %   blk = env.AgentBlock
        AgentBlock string
    end
    properties
        % ENV.RESETFCN
        %   Specify a function handle that will be executed before each
        %   simulation. The function handle must take in and return a
        %   scalar Simulink.SimulationInput object.
        %
        %   env.ResetFcn = @(in) setVariable(in,'x0',rand(1,1));
        %
        % See also SIMULINK.SIMULATIONINPUT
        ResetFcn = []
        
        % ENV.USEFASTRESTART
        %   Toggle wheter to use fast restart mode when simulating/training
        %   the environment model (default is 'on')
        %
        %   env.UseFastRestart = 'on';
        UseFastRestart = 'on'
    end
    properties (Access = private,Transient)
        
        % temporary/user set policy variable
        TemporaryPolicyVariable
        OriginalPolicyVariable
        
        % simulink simulation engine
        SimEngine = []
        SimMgr    = [];
    end
    methods
        %% set/get methods
        function set.UseFastRestart(this,val)
            % set the fast restart property
            validatestring(val,{'on','off'},'','UseFastRestart');
            this.UseFastRestart = val;
        end
        function set.ResetFcn(this,fcn)
            % install the reset function
            if ~isempty(fcn)
                validateattributes(fcn,{'function_handle'},{'scalar'},'','ResetFcn');
                if nargin(fcn) ~= 1
                    error(message('rl:env:SimulinkBlockEnvInvalidResetFcn'));
                end
            end
            this.ResetFcn = fcn;
        end
        %% constructor
        function this = SimulinkEnvWithAgent(mdl,blk,observationInfo,actionInfo,varargin)
            % parse arguments
            p = inputParser;
            addRequired(p,'model'     ,@(x) validateattributes(x,{'char','string'},{'scalartext'}));
            addRequired(p,'agentBlock',@(x) validateattributes(x,{'char','string'},{'scalartext'}));
            % error checking for obs/act info will happen inside the env
            % abstraction
            addRequired(p,'observationInfo',@(x) validateattributes(x,{'rl.util.RLDataSpec'},{},'','obervationInfo'));
            addRequired(p,'actionInfo'     ,@(x) validateattributes(x,{'rl.util.RLDataSpec'},{},'','actionInfo'    ));
            addParameter(p,'UseFastRestart','on');
            % parse
            parse(p,mdl,blk,observationInfo,actionInfo,varargin{:});
            % extract
            mdl = p.Results.model;
            blk = p.Results.agentBlock;
            observationInfo = p.Results.observationInfo;
            actionInfo = p.Results.actionInfo;
            usefastrestart = p.Results.UseFastRestart;
            
            this = this@rl.env.AbstractEnv(observationInfo,actionInfo);
            
            % make sure the model is load-able
            allmdls = localLoadModel(mdl);
            
            % make sure the root model is consistent with the block
            % path
            allagents = {};
            for i = 1:numel(allmdls)
                % make sure the referenced models are loaded
                % make sure the model is load-able
                allagents = union(allagents,find_system(allmdls{i},'LookUnderMasks','on','FollowLinks','on','MaskType','RL Agent'));
            end
            
            % make sure the agent block is our rl library block
            if ~strcmp(get_param(blk,'MaskType'),'RL Agent')
                error(message('rl:env:SimulinkBlockEnvInvalidAgentBlock'));
            end
            
            % make sure the agent block is a member of the root model
            if ~ismember(blk,allagents)
                error(message('rl:env:SimulinkBlockEnvAgentBlockNotMemberOfModel',blk,mdl));
            end
            
            this.Model = string(mdl);
            this.AgentBlock = string(blk);
            
            % use fast restart
            this.UseFastRestart = usefastrestart;
        end
        %% destructor
        function delete(~)
            % destroy the environment
            
            % Nothing to do here as all properties are managed at sim time
            % (see simWrapper and parsimWrapper)
        end
        %% methods needed by the agent
        function [nextobs,reward,isdone,loggedSignals] = step(this,actions) %#ok<STOUT,INUSD>
            % step the environment given an action
            
            % step is not yet supported for simulink
            error(message('rl:env:SimulinkBlockEnvStepNotSupported'));
        end
        
        function in = reset(this)
            % IN = RESET(ENV)
            %   Return a Simulink.SimulationInput object based ResetFcn if
            %   provided. If ResetFcn is empty, a default
            %   Simulink.SimulationInput object will be returned.
            %
            % See also SIMULINK.SIMULATIONINPUT
            in = getSimulationInput(this);
            in = localResetOnSimInput(this.ResetFcn,in);
        end
        function stop(this)
            % STOP(ENV)
            %   stop any running simulations on the environment or
            %   terminate the model if it is being stepped.
            if isModelLoaded(this)
                set_param(this.Model,'SimulationCommand','stop');
            end
        end
        function validateEnvironment(this)
            % Validate Simulink environment by performing a short simulation for 2 sample times.
            % The validation assumes the agent is already created in the base workspace.
            
            % collect policy info
            blk = this.AgentBlock;
            mdl = this.Model;
            varname = get_param(blk,'Agent');
            if Simulink.data.existsInGlobal(mdl,varname)
                currentPolicy = Simulink.data.evalinGlobal(mdl,varname);
            else
                error(message('rl:env:errSimulinkEnvCheckNoAgent'));
            end
            % ensure agent only sim during validation: set step mode to "sim"
            currentStepMode = getStepMode(currentPolicy);
            setStepMode(currentPolicy,"sim");       
            clnup = onCleanup(@() setStepMode(currentPolicy,currentStepMode));
            
            % configure to run sim up to 2 time steps
            Ts = getSampleTime(currentPolicy);
            T0 = Simulink.data.evalinGlobal(mdl,get_param(mdl,'StartTime'));
            Te = T0 + 2*Ts; % take 2 steps
            
            % create the sim function
            resetfcn = this.ResetFcn;
            
            simData.Model = mdl;
            simData.Block = blk;
            simData.ResetFcn = resetfcn;
            simData.EpisodeTime = Te;
            simData.InitialTime = T0;
            simData.UsePCT = false;

            % validate step function, utilize try/catch in step()
            simout = localInnerSimFcn(simData,getSimulationInput(this),1); %#ok<NASGU>
        end
    end
    methods (Access = private)
        function t = getSimTime(this)
            % get the time at the current sample hit
            errorIfNotLoaded(this,'getSimTime');
            t = get_param(this.Model,'simulationtime');
        end
        function errorIfNotLoaded(this,caller)
            if ~isModelLoaded(this)
                error(message('rl:env:SimulinkBlockEnvModelMustBeLoadedToCall',this.Model,caller));
            end
        end
        function val = isModelLoaded(this)
            % check if the model is loaded
            val = ~isempty(this.Model) && bdIsLoaded(this.Model);
        end
        function prepare(this,policy)
            
            % resolve the policy
            resolvePolicy(this,policy);
            
            mdl = this.Model;
            
            % the rl agent block is only supported in normal and
            % accelerator mode
            mode = get_param(mdl,'SimulationMode');
            if ~ismember(mode,{'normal','accelerator'})
                error(message('rl:env:SimulinkBlockEnvInvalidSimMode'));
            end
        end
        function restore(this)
            % restore the model to it's original settings
                        
            % make sure we stop any type of simulation
            stop(this);
            
            % restore the policy
            restorePolicy(this);
        end
        function in = getSimulationInput(this)
            % convenience function for generating a Simulink simulation
            % input
            in = Simulink.SimulationInput(this.Model);
        end
        function cancelSims(this)
            if ~isempty(this.SimEngine)
                cancel(this.SimEngine);
            end
        end
        function cancelSimsByPolicy(this,policy)
            if policy.TerminateSimulation
                cancelSims(this);
            end
        end
        function simouts = executeSimsWrapper(this,policy,in,simfh,simouts,opts)
            % wrapper around executeSims on the simulink engine
            
            
            %% nested fcn for gathering sim outputs
            finishedSimCount = 0;
            ex_ = [];
            function nestedSimFinishedBC(ed)
                % get the simoutput
                simout_ed = ed.SimulationOutput;
                % g2153021
                if isprop(simout_ed,'RL_TEMP_EpisodeInfo')
                    finishedSimCount = finishedSimCount + 1;
                    % handle any sim errors
                    try
                        localHandleSimoutErrors(simout_ed,opts.StopOnError);
                    catch ex_
                        % cancel the remaining sims
                        cancelSims(this);
                    end
                    % tell the world an episode has finished
                    if ~policy.TerminateSimulation
                        simouts{finishedSimCount} = simout_ed;
                        notifyEpisodeFinished(this,...
                            simout_ed.RL_TEMP_EpisodeInfo,...
                            simout_ed,...
                            finishedSimCount,...
                            simout_ed.RL_TEMP_WorkerID);
                    end
                end
            end
            
            %% setup listeners on the engine
            simlist(1) = event.listener(this.SimMgr,'SimulationFinished'   ,@(src,ed) nestedSimFinishedBC(ed));
            simlist(2) = event.listener(this.SimMgr,'AllSimulationsQueued' ,@(src,ed) notifySimsDistributed(this));
            simlist(3) = addlistener(policy,'TerminateSimulation','PostSet',@(src,ed) cancelSimsByPolicy(this,policy));
            
            clnup = onCleanup(@() delete(simlist));
            
            % turn off executeSims related warnings if the errors will be
            % thrown anyway
            ws = warning('off','Simulink:Commands:SimulationsWithErrors');
            wcleanup = onCleanup(@() warning(ws));
            
            %% simulate
            try
                executeSims(this.SimEngine,simfh,in);
            catch ex
                % ignore abort sim errors
                ids2ignore = ["Simulink:Commands:SimAborted","parallel:fevalqueue:ExecutionCancelled"];
                if ~ismember(ex.identifier,ids2ignore)
                    rethrow(ex);
                end
            end
            %% throw errors
            localReframeSimError(this.Model,this.AgentBlock,ex_);
        end
        function simout = simWrapper(this,policy,simData,in,opts)
            % wrapper of sim with Simulink.SimulationInput object
            
            % resolve the policy
            prepare(this,policy);
            clnup = onCleanup(@() restore(this));
            
            N = opts.NumSimulations;
            simouts = cell(N,1);
                       
            simfh = @(in_) localInnerSimFcn(simData,in_,[]);
            simouts = executeSimsWrapper(this,policy,in,simfh,simouts,opts);
            
            simout = vertcat(simouts{:});
        end
        function simout = parsimWrapper(this,policy,simData,in,opts)
            % wrapper of parsim with Simulink.SimulationInput object
            
            % resolve the policy
            cenv = parallel.pool.Constant(this);
            localPrepareOnWorker(cenv,policy);
            clnup = onCleanup(@()localRestoreOnWorker(cenv));
            
            N = opts.NumSimulations;
            simouts = cell(N,1);
                
            % create a constant for the simData
            csimData = parallel.pool.Constant(simData);

            % create the fh
            simfh = @(in_) localInnerSimFcn(csimData.Value,in_,[]);

            % run in background. executeSims will give us the future
            if this.SimEngine.Options.RunInBackground

                % call parfeval (via executeSims)
                F = executeSims(this.SimEngine,simfh,in);

                % tell the world that the sims have been distributed
                notifySimsDistributed(this);

                for simCount = 1:N
                    if ~policy.TerminateSimulation
                        [~,simout_] = fetchNext(F);
                        
                        try
                            localHandleSimoutErrors(simout_,opts.StopOnError);
                        catch ex
                            localReframeSimError(this.Model,this.AgentBlock,ex);
                        end
                        
                        simouts{simCount} = simout_;
                        % tell the world an episode has finished
                        notifyEpisodeFinished(this,...
                            simout_.RL_TEMP_EpisodeInfo,simout_,simCount,simout_.RL_TEMP_WorkerID);
                    else
                        cancel(F);
                        break;
                    end
                end
            else
                simouts = executeSimsWrapper(this,policy,in,simfh,simouts,opts);
            end
            simout = vertcat(simouts{:});
        end
        function simCleanup(this)
            % restore the model to it's original state
            restore(this);
        end
        function resolvePolicy(this,policy)
            % make sure the input policy is the same policy used for
            % simulation in simulink. It it's not, create a temporary
            % variable in the global workspace
            
            blk = this.AgentBlock;
            mdl = this.Model;
            
            % localGetModelPolicy will load the model
            [currentPolicy,varname] = localGetModelPolicy(mdl,blk);
                    
            if ~isequal(policy,currentPolicy)
                % check if model is already compiled, as we have to turn
                % fast restart off to set the policy on the block
                if strcmp(get_param(mdl,'FastRestart'),'on')
                    % REVISIT throw a warning?
                    set_param(mdl,'FastRestart','off');
                    clnup = onCleanup(@() set_param(mdl,'FastRestart','on'));
                end
                
                % push the policy into a static ws
                rl.util.PolicyInstance.set(policy);
                tempname = 'rl.util.PolicyInstance.get()';
                
                this.TemporaryPolicyVariable = tempname;
                this.OriginalPolicyVariable  = varname ;
                set_param(blk,'Agent',tempname);
            end
        end
        function restorePolicy(this)
            % restore the user set policy on the block
            varname = this.OriginalPolicyVariable;
            mdl = this.Model;
            if ~isempty(varname)   
                wasloaded = bdIsLoaded(mdl);
                if wasloaded
                    % temporarily turn FR off so we can set the agent param
                    if strcmp(get_param(mdl,'FastRestart'),'on')
                        % REVISIT throw a warning?
                        set_param(mdl,'FastRestart','off');
                        clnup = onCleanup(@() set_param(mdl,'FastRestart','on'));
                    end
                    blk = this.AgentBlock;
                    set_param(blk,'Agent',varname);
                end
                
                % flush the static ws
                rl.util.PolicyInstance.set([]);

                % close the model again if it wasn't loaded
                if ~wasloaded
                    close_system(mdl,0);
                end
            end
            this.TemporaryPolicyVariable = [];
            this.OriginalPolicyVariable  = [];
        end
    end
    methods (Hidden)
        function name = getNameForEpisodeManager(this)
            name = this.Model;
        end
    end
    methods (Access = protected)
        function [experiences,varargout] = simWithPolicyImpl(env,policy,opts,varargin)
            % simulate the environment/Simulink model. This defines the
            % implementation of simWithPolicy
                        
            % if an agent and max steps is provided, overload final
            % simulation time. Note, Simulink will call getInitialAction at
            % the first sample hit, then will call step on all subsequent
            % sample hits
            Ts = getSampleTime(policy);
            Te = Ts*(opts.MaxSteps);
            T0 = Simulink.data.evalinGlobal(env.Model,get_param(env.Model,'StartTime'));
            
            % build the simulation input objects
            N = opts.NumSimulations;
            in(1:N) = getSimulationInput(env);
            
            % create the sim function
            mdl = env.Model;
            blk = env.AgentBlock;
            resetfcn = env.ResetFcn;
            usePCT = opts.UseParallel && matlab.internal.parallel.isPCTInstalled();
            
            simData.Model = mdl;
            simData.Block = blk;
            simData.ResetFcn = resetfcn;
            simData.EpisodeTime = Te;
            simData.InitialTime = T0;
            simData.UsePCT = usePCT;
            
            try
                % don't check worker id here g2036797
                if usePCT
                    % parsim
                    simouts = parsimWrapper(env,policy,simData,in,opts);
                else
                    % sim
                    simouts = simWrapper(env,policy,simData,in,opts);
                end
            catch ex
                throwAsCaller(ex);
            end
            
            if isempty(simouts)
                experiences = {};
            else
                experiences = {simouts.RL_TEMP_ExperienceCell}';
                % remove RL_ExperienceCell and RL_EpisodeInfo from the output
                simouts = removeProperty(simouts,'RL_TEMP_ExperienceCell'); 
                simouts = removeProperty(simouts,'RL_TEMP_EpisodeInfo'   ); 
                simouts = removeProperty(simouts,'RL_TEMP_WorkerID'      ); 
            end
            
            if nargout > 1
                varargout{1} = simouts;
            end
        end
        function setupForSimImpl(env,opts)
            
            % make sure the sim engine is initialized
            % generate a "dummy" SimulationInput object
            % in_(1:opts.NumSimulations) = reset(env);
            if isempty(env.SimEngine)
                in_ = getSimulationInput(env);

                % get the sim manager
                env.SimMgr = Simulink.SimulationManager(in_);

                % get the sim engine
                env.SimEngine = env.SimMgr.SimulationManagerEngine;
            end
            simEng = env.SimEngine;

            % common options between parallel and local sim
            simEng.Options.UseFastRestart = strcmpi(env.UseFastRestart,"on");
            % Always store the error as part of the simoutput. We will
            % throw the error if StopOnError is "on" once simoutput is
            % received
            simEng.Options.StopOnError                      = false;
            simEng.Options.ShowSimulationManager            = false;
            simEng.Options.ShowProgress                     = false;

            if opts.UseParallel
                paropts    = opts.ParallelizationOptions;
                afiles     = paropts.AttachedFiles;
                setupFcn   = paropts.SetupFcn;
                cleanupFcn = paropts.CleanupFcn;
                xferbasews = paropts.TransferBaseWorkspaceVariables;

                if isempty(afiles)
                    afiles = {};
                else
                    afiles = cellstr(afiles);
                end
                % determine the necessary files to run the reset function
                resetfiles = localGetFilesFromFunctionHandles({env.ResetFcn});
                if isempty(afiles)
                    afiles = resetfiles;
                else
                    afiles = vertcat(afiles(:),resetfiles(:));
                end

                % configure the engine options
                simEng.Options.SetupFcn                         = setupFcn;
                simEng.Options.CleanupFcn                       = cleanupFcn;
                simEng.Options.UseParallel                      = true;
                % run in background is not needed since the simMgr
                % has the events AllSimulationsQueued and
                % SimulationFinished
                simEng.Options.RunInBackground                  = false;
                simEng.Options.AttachedFiles                    = afiles;
                simEng.Options.TransferBaseWorkspaceVariables   = strcmpi(xferbasews,"on");
            else
                simEng.Options.SetupFcn                         = [];
                simEng.Options.CleanupFcn                       = [];
                simEng.Options.UseParallel                      = false;
                simEng.Options.RunInBackground                  = false;
                simEng.Options.AttachedFiles                    = {};
                simEng.Options.TransferBaseWorkspaceVariables   = false;
            end

            % run the setup
            setup(simEng);
        end
        function cleanupForSimImpl(env)              
            % cleanup from the sim engine
            if ~isempty(env.SimEngine)
                cleanup(env.SimEngine);
            end
            delete(env.SimEngine);
            delete(env.SimMgr);
            env.SimEngine = [];
            env.SimMgr    = [];
        end
        function policy = getAttachedPolicyImpl(env)
            policy = localGetModelPolicy(env.Model,env.AgentBlock);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Local Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function simout = localInnerSimFcn(simData,in,simCount)

% this function allows us to "lazily" reset and setup the sim
% inputs which can be expensive for lots of episodes

mdl           = simData.Model;
blk           = simData.Block;
resetfcn      = simData.ResetFcn;
Te            = simData.EpisodeTime;
T0            = simData.InitialTime;
usePCT        = simData.UsePCT;

if isempty(simCount)
    simCount = in.RunId;
end

% supports parallel
policy = localGetModelPolicy(mdl,blk);

% call reset on the sim input
in = localResetOnSimInput(resetfcn,in);

% get the final time, if there is an operating point use that as T0,
% otherwise use the model StartTime. Note, this has to happen after the
% reset function in the case the user installs a State
if isempty(in.InitialState)
    Tf = T0 + Te;
else
    Tf = in.InitialState.snapshotTime + Te;
end
Tf_str = num2str(Tf);

% set the final time
in = setModelParameter(in,'StopTime',Tf_str);

% capture errors
in = addHiddenModelParameter(in,'CaptureErrors','on');

% execute the policy pre sim fcn
preSimFcn(policy,simCount);

% call simulink sim
simout = sim(in);

% execute the policy post sim fcn
postSimFcn(policy,simCount,simout);

% attach the experiences and episode info to the simulation
% output
simout.RL_TEMP_ExperienceCell = getExperiences(policy);
simout.RL_TEMP_EpisodeInfo    = getEpisodeInfo(policy);
simout.RL_TEMP_WorkerID       = rl.util.getWorkerID(usePCT);
end
function [currentPolicy,varname] = localGetModelPolicy(mdl,blk)
% get the policy configured in the model

% make sure the model is loaded before get_param is called
% g1967773
localLoadModel(mdl);

% get the var name from the block
varname = get_param(blk,'Agent');

% get the policy from the global ws
try
    currentPolicy = Simulink.data.evalinGlobal(mdl,varname);
catch
    currentPolicy = [];
end
end
function allMdls = localLoadModel(mdl)
% load the model and throw any errors

try
    load_system(mdl);
    allMdls = find_mdlrefs(mdl);
catch ex
    me = MException(message('rl:env:SimulinkBlockEnvCouldNotLoadModel',mdl));
    me = addCause(me,ex);
    throw(me);
end
% load referenced models
for ct = 1:numel(allMdls)
    rmdl = allMdls{ct};
    if ~strcmp(mdl,rmdl)
        try
            load_system(rmdl);
        catch ex
            me = MException(message('rl:env:SimulinkBlockEnvCouldNotLoadRefModel',allMdls{ct}));
            me = addCause(me,ex);
            throw(me);
        end
    end
end
end
function in = localResetOnSimInput(resetfcn,in)
% call reset on the sim input
if ~isempty(resetfcn)
    try
        in = resetfcn(in);
        validateattributes(in,{'Simulink.SimulationInput'},{'scalar'});
    catch ex
        me = MException(message('rl:env:SimulinkBlockEnvErrorEvaluatingResetFcn'));
        me = addCause(me,ex);
        throw(me);
    end
end
end
function files = localGetFilesFromFunctionHandles(functionHandles)
% get the files needed to execute the code specified by function handles
functionFiles = MultiSim.internal.getFilenamesForFunctionHandles(functionHandles);
functionFiles = functionFiles(~strcmp(functionFiles,''));
files = matlab.codetools.requiredFilesAndProducts(functionFiles);
end
function localHandleSimoutErrors(simout,stopOnError)
if strcmp(stopOnError,"on")
    ediag = simout.SimulationMetadata.ExecutionInfo.ErrorDiagnostic;
    if ~isempty(ediag)
        ex = ediag.Diagnostic;
        ids2ignore = ["Simulink:Commands:SimAborted","parallel:fevalqueue:ExecutionCancelled"];
        if ~ismember(ex.identifier,ids2ignore)
            reportAsError(ex);
        end
    end
end
end
function localReframeSimError(mdl,blk,ex)
if ~isempty(ex)
    % generic error for multiple cases
    me = MException(message('rl:env:SimulinkBlockEnvErrorDuringSimulation',...
        mdl,get_param(blk,'Agent')));
    switch ex.identifier
        case {...
                'SystemBlock:MATLABSystem:MethodInvokeError',...
                'SystemBlock:MATLABSystem:MethodInvokeErrorWithStack'...
                }
            me = addCause(me,ex.cause{1});
        case 'Simulink:Parameters:NoChangeWhenInFastRestart'
            me = MException(message('rl:env:SimulinkBlockEnvErrorWithFastRestartOn'));
            me = addCause(me,ex);
        otherwise
            me = addCause(me,ex);
    end
    throw(me);
end
end
function localPrepareOnWorker(this_const,policy)
% will make sure the policy is built on the workers
F = parfevalOnAll(@() prepare(this_const.Value,policy),0);
% handle any errors
try
    fetchOutputs(F);
catch ex
    throwAsCaller(ex);
end
end
function localRestoreOnWorker(this_const)
% will restore the policy on the worker
F = parfevalOnAll(@() restore(this_const.Value),0);
% handle any errors
try
    fetchOutputs(F);
catch ex
    throwAsCaller(ex);
end
end