classdef (Abstract) AbstractEnv < handle
    % RLENVABSTRACT: Creates abstract class for Reinforcement Learning
    % environment.

    % Copyright 2017-2018 The MathWorks Inc.
    
    properties (Access = protected)
        % Information on observations
        ObservationInfo rl.util.RLDataSpec
        % Information on actions
        ActionInfo rl.util.RLDataSpec
    end
    properties (Access = protected,Transient)
        % options cached by setupForSim. Cleaned up by cleanupForSim.
        SetupForSimOptions = [];
    end
    events
        % tell the world an episode has been finished
        EpisodeFinished
        % tell the world the set of simulations has been distributed
        SimulationsDistributed
    end
    methods
        function this = AbstractEnv(observationInfo,actionInfo)
            % Abstract environment construction, requires observation and
            % action info
            try 
                narginchk(2,2)
            catch
                error(message('rl:env:ActionAndObservationInfoNotDefined'));
            end
            
            this.ObservationInfo = observationInfo;
            this.ActionInfo      = actionInfo     ;       
        end
    end
    
    methods (Abstract)
        
        % We may want to use LoggedStates for calulating reward (outside
        % the environment) on signals that are not necessarily provided as
        % observations
        
        % [OBSERVATION,REWARD,ISDONE,LOGGEDSIGNALS] = step(ENV,ACTION) 
        % simulates the environment ENV with given ACTION for one step.
        % step returns OBSERVATION, REWARD, ISDONE terminal flag and any
        % additional information LOGGEDSIGNALS.
        [obsevation,reward,isTerminal,loggedSignals] = step(this,action)
        
        % OBSERVATION = reset(ENV) returns the initial OBSERVATION of the ENV.
        InitialState = reset(this)
    end
    methods (Hidden)
        % simulate the environment against a policy give a set of
        % simulation options
        function [experiences,varargout] = simWithPolicy(this,policy,opts,varargin)
            % push the max steps to the policy
            policy.MaxSteps = opts.MaxSteps;
            
            % setup for sim (NOTE: the env can be setup external to
            % simWithPolicy however it is up to the user to cleanup
            % properly)
            if ~isSetupForSim(this)
                setupForSim(this,opts);
                % cleanup the workers at the end
                clnup = onCleanup(@() cleanupForSim(this));
            end
            
            % run the subclass method
            [experiences,varargout{1:(nargout-1)}] = simWithPolicyImpl(this,policy,opts,varargin{:});
        end
    end
    methods (Access = protected,Abstract)
        [experiences,varargout] = simWithPolicyImpl(this,policy,opts,varargin)
    end
    methods        
        function Info = getObservationInfo(this)
            Info = this.ObservationInfo;
        end        
        function Info = getActionInfo(this)
            Info = this.ActionInfo;         
        end
        % Simulates environment and agent
        function varargout = sim(env,policy,varargin)
            % EXPERIENCES = SIM(ENVIRONMENT,AGENT,OPTIONS) 
            % EXPERIENCES = SIM(AGENT,ENVIRONMENT,OPTIONS)
            %
            % Simulate the ENVIRONMENT against an AGENT with simulation options OPTIONS
            % created from rlSimulationOptions.  
            %
            % An EXPERIENCES structure is returned with the number of elements equal to
            % the number of simulations run. The EXPERIENCES structure has the
            % following fields:
            %
            %   Observation         : Structure with fields corresponding to the
            %                         observation info of the environment. Each field
            %                         is a timeseries of the respective observation
            %                         signal.
            %   Action              : Structure with fields corresponding to the
            %                         action info of the environment. Each field
            %                         is a timeseries of the respective action
            %                         signal.
            %   Reward              : timeseries of the scalar reward signal.
            %   IsDone              : timeseries of the scalar IsDone signal,
            %                         signifying the end of episode.
            %
            % Example: Simulate the Double Integrator environment against an agent for
            % up to 1000 steps.
            %   
            %   env = rlPredefinedEnv('DoubleIntegrator-Discrete'); 
            %   opt = rlSimulationOptions('MaxSteps',1000); 
            %   experiences = sim(env,agent,opt);
            %
            % See also: rlSimulationOptions, timeseries
            
            % parse the input arguments
            p = inputParser;
            addRequired(p,'environment',@(x)validateattributes(x,{'rl.env.AbstractEnv'},{'scalar'}));
            addRequired(p,'policy',@(x)validateattributes(x,{'rl.agent.AbstractPolicy'},{'scalar'}));
            addOptional(p,'options',rlSimulationOptions(),@(x)validateattributes(x,{'rl.option.rlSimulationOptions'},{'scalar'}));
            
            parse(p,env,policy,varargin{:});
            opts = p.Results.options;
            
            % collect experiences if outputs are requested
            collectExperiences = nargout > 0;
            if collectExperiences
                attachLogger(policy,opts.MaxSteps);
            end
            
            % make sure the policy is using sim mode
            setStepMode(policy,"sim");
            
            % simulate
            [experiencesCell,simInfo] = simWithPolicy(env,policy,opts);
            
            % create the output experiences if requested
            if isempty(experiencesCell) || isempty(experiencesCell{1})
                % if all errors occur before elements can be logged,
                % created an empty structure with the sim info which will
                % contain the error.
                for i = numel(simInfo):-1:1
                    s(i).Observation     = [];
                    s(i).Action          = [];
                    s(i).Reward          = [];
                    s(i).NextObservation = [];
                    s(i).IsDone          = [];
                    s(i).SimulationInfo  = simInfo(i);
                end
                experiences = s(:);
            else
                % get the sample time from the policy
                Ts = getSampleTime(policy);
                % convert Ts of -1 to 1
                if Ts == -1
                    Ts = 1;
                end
                N = size(experiencesCell,1);
                % store the individual sims in a cell for easy
                % concatenation later on
                outputExperienceCell = cell(N,1);
                ainfo = getActionInfo     (env);
                oinfo = getObservationInfo(env);
                
                odatatype = cellstr([oinfo.DataType]);
                adatatype = cellstr([ainfo.DataType]);
                for i = 1:N
                    experience = experiencesCell{i};
                    
                    observation     = cellfun(@(x)x{1},experience,'UniformOutput',false);
                    action          = cellfun(@(x)x{2},experience,'UniformOutput',false);
                    reward          = cellfun(@(x)x{3},experience,'UniformOutput',false);
                    nextobservation = cellfun(@(x)x{4},experience,'UniformOutput',false);
                    isdone          = cellfun(@(x)x{5},experience,'UniformOutput',false);
                    
                    % push the last observation from nextobservation onto
                    % observation
                    obs_ = [observation;nextobservation(end)];
                    
                    for j = size(obs_{1},2):-1:1
                        obs{j} = cellfun(@(x) x{j},obs_,'UniformOutput',false);
                        
                        % get the dimension of the observation
                        d = obs{j}{1};
                        nd = ndims(d);
                        
                        % cat along the last dimension + 1
                        obs{j} = cat(nd+1,obs{j}{:});
                        
                        % cast the observation to the appropriate data type
                        if isnumeric(obs{j})
                            obs{j}  = cast(obs{j},odatatype{j});
                        end
                    end
                    for j = size(action{1},2):-1:1
                        act{j} = cellfun(@(x) x{j},action,'UniformOutput',false);
                        
                        % get the dimension of the action
                        d = act{j}{1};
                        nd = ndims(d);
                        
                        % cat along the last dimension + 1
                        act{j} = cat(nd+1,act{j}{:});
                        
                        % cast the action to the appropriate data type
                        if isnumeric(act{j})
                            act{j} = cast(act{j},adatatype{j});
                        end
                    end
                    
                    rwd_ = [reward{:}]';
                    isd_ = [isdone{:}]';
                    
                    % set non-terminal conditions back to zero
%                     isd_(isd_ > 1) = 0;
                        
                    timevec = (0:(numel(rwd_)-1))*Ts;
                    
                    % create reward and isdone timeseries
                    rwd  = timeseries(rwd_,timevec,'Name','Reward');
                    isd  = timeseries(isd_,timevec,'Name','IsDone');
                    
                    % create the timeseries objects from the IO data
                    for j = 1:numel(oinfo)
                        name = oinfo(j).Name;
                        if isempty(name) || name == ""
                            name = sprintf('obs%u',j);
                        end
                        vname = matlab.lang.makeValidName(name);
                        sobs.(vname) = timeseries(obs{j},[timevec,timevec(end)+Ts],'Name',name);
                    end
                    for j = 1:numel(ainfo)
                        name = ainfo(j).Name;
                        if isempty(name) || name == ""
                            name = sprintf('act%u',j);
                        end
                        vname = matlab.lang.makeValidName(name);
                        sact.(vname) = timeseries(act{j},timevec,'Name',name);
                    end
                    
                    % store in an output structure
                    s.Observation     = sobs;
                    s.Action          = sact;
                    s.Reward          = rwd;
                    s.IsDone          = isd;
                    s.SimulationInfo  = simInfo(i);
                    
                    outputExperienceCell{i} = s;
                end
                % concat experiences
                experiences = vertcat(outputExperienceCell{:});
            end
            
            if collectExperiences
                % only return experiences if requested, since if a logger
                % is not attached, all experiences will be empty
                varargout{1} = experiences;
            end
            
            % detatch logger
            detatchLogger(policy);
        end
        % Validate the environment by performing a short simulation
        function validateEnvironment(this) %#ok<MANU>
        end
    end
    methods (Abstract,Access = protected)
        setupForSimImpl(this,opts)
        cleanupForSimImpl(this)
    end
    methods (Hidden)
        function name = getNameForEpisodeManager(this)
            name = regexprep(class(this),'\w*\.','');
        end
        function setupForSim(this,opts)  
            % setup the env for simulation
            
            % if already setup, cleanup first
            if isSetupForSim(this)
                cleanupForSim(this);
            end
            
            % setup the subclass env
            setupForSimImpl(this,opts);
            
            % if running in parallel, set the random seed
            if opts.UseParallel
                wseeds = opts.ParallelizationOptions.WorkerRandomSeeds;
                pool = gcp();
                numWorkers = pool.NumWorkers;
                
                % get the ids from the workers g2036797
                F = parfevalOnAll(@rl.util.getWorkerID,1);
                ids = fetchOutputs(F);
                
                % expand the worker random seeds
                if isscalar(wseeds)
                    wseeds = wseeds*ones(1,numWorkers);
                end
                if numel(wseeds) ~= numWorkers
                    error(message('rl:general:ParallelTrainInvalidNumberOfWorkerSeeds',numel(wseeds),numWorkers));
                end
                
                % assign the seeds to the worker ids
                seedmap = containers.Map('KeyType','int32','ValueType','int32');
                for i = 1:numWorkers
                    seedmap(ids(i)) = wseeds(i);
                end
                
                % setup the seeds
                F = parfevalOnAll(@localSetupRandomSeedOnWorkers,0,seedmap);
                
                % handle any errors
                try
                    fetchOutputs(F);
                catch ex
                    throwAsCaller(ex);
                end
            end
            
            % cache the options
            this.SetupForSimOptions = opts;
        end
        function cleanupForSim(this)
            if isSetupForSim(this)
                cleanupForSimImpl(this)
                % cleanup the options
                this.SetupForSimOptions = [];
            end
        end
        function val = isSetupForSim(this)
            % check if the env is setup for sim
            val = ~isempty(this.SetupForSimOptions);
        end
    end
    methods (Access = protected)
        function notifySimsDistributed(this)
            % tell the world that simulations have been distributed
            notify(this,'SimulationsDistributed');
            drawnow();
        end
        function notifyEpisodeFinished(this,episodeInfo,simInfo,simCount,workerID)
            % tell the world that an episode has been finished and provide
            % the episode info/sim #
            if nargin < 5
                workerID = 0;
            end
            info.EpisodeInfo    = episodeInfo;
            info.EpisodeCount   = simCount;
            info.WorkerID       = workerID;
            info.SimulationInfo = simInfo;
            ed = rl.util.RLEventData(info);
            notify(this,'EpisodeFinished',ed);
            drawnow();
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Local Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localSetupRandomSeedOnWorkers(seedmap)
id = rl.util.getWorkerID();
ws = seedmap(id);
if ws ~= -2
    % if set to default (-1) use the id for the seed
    if ws == -1
        ws = id;
    end
    rng(ws);
end
end