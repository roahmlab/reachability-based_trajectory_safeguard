classdef ParforTrainTask < rl.task.Task
% PARFORTRAINTASK
%
% Run training using parfor. Useful for sync parallel methods that can wait
% to train until the end of an episode.

% Revised: 8-14-2019
% Copyright 2019 The MathWorks, Inc.

    properties
        TrainOpts
        SimOpts
        % the agent that ends up being distributed to the workers
        DistributedAgent
        Env
    end
    properties (Transient)
        Agent
    end
    events
        WorkersSetup
        WorkersCleanedup
        EpisodeInfoReceivedFromWorker
        ActionMessageReceived
    end
    methods
        function this = ParforTrainTask(env,agent,trainOpts)
            this = this@rl.task.Task();
            
            this.Agent = agent;
            this.Env = env;
            this.TrainOpts = trainOpts;
            
            simOpts = getSimulationOptions(trainOpts);
            simOpts.UseParallel = false;
            simOpts.NumSimulations = 1;
            
            this.SimOpts = simOpts;
            
            % create a "distributed copy" of the agent
            this.DistributedAgent = buildDistributedAgent(this);
        end
    end
    methods
        function [workerData,simInfo,workerEpisodeInfo] = parallelExecution(this,params,isFirstLoop)
            % code that is run inside parfor
            
            if isFirstLoop
                setupParforWorkers(this);
            end
            
            % set the updated tunable parameters
            setLearnableParameters(this.DistributedAgent,params);
            
            % run the sim and capture data for learning
            [workerData,simInfo,workerEpisodeInfo] = runWorkerSims(this);
        end
        function learnFromWorkerData(this,wdata,numSims)
            learnFromWorkerDataImpl(this,wdata,numSims);
        end
        function str = generateDataReceivedMsg(this,wdata,wid)
            str = generateDataReceivedMsgImpl(this,wdata,wid);
        end
    end
    methods (Abstract,Access = protected)
        setupParforWorkers(this)
        [wdata,wsimInfo,wepInfo] = runWorkerSims(this)
        learnFromWorkerDataImpl(this,wdata,numSims)
        str = generateDataReceivedMsgImpl(this,wdata,wid)
    end
    methods (Access = protected)
        
        function agent = buildDistributedAgent(this)
            agent = this.Agent;
        end
        
        function varargout = runImpl(this)
            varargout = {};
            
            % setup the pool if not already setup
            pool = gcp();
            numWorkers = pool.NumWorkers;
            
            trainOpts = this.TrainOpts;
            simOpts = this.SimOpts;
            env = this.Env;
            agent = this.Agent;
            
            % keep training until there are no more remaining sims
            remainingSims = trainOpts.MaxEpisodes;
            
            % distribute the env across workers
            setupForSim(env,simOpts);
            clnup1 = onCleanup(@()cleanupForSim(env));
            clnup2 = onCleanup(@()setStepMode(agent,"sim"));
            
            % notify sims distributed
            notify(this,"WorkersSetup");
            clnup3 = onCleanup(@()notify(this,"WorkersCleanedup"));
            
            % create worker constants to avoid copies at the start of each
            % parfor
            cthis = parallel.pool.Constant(this);
            
            isFirstLoop = true;
            episodeCount = 0;
            while remainingSims > 0 && ~agent.TerminateSimulation
                
                % figure out how many sims to run (up to the number of available
                % workers)
                activeSims = min(remainingSims,numWorkers);
                % create storage vars
                workerData        = cell(1,activeSims);
                workerEpisodeInfo = cell(1,activeSims);
                simInfos          = cell(1,activeSims);
                
                % get the updated tunable params
                params = getLearnableParameters(agent);
                
                % run the simulations on the workers
                parfor i = 1:activeSims
                    [workerData{i},simInfos{i},workerEpisodeInfo{i}] = parallelExecution(cthis.Value,params,isFirstLoop)
                end
                % learn from the gathered data
                for i = 1:activeSims
                    if agent.TerminateSimulation
                        break;
                    end
                    
                    % notify data received
                    ed = rl.util.RLEventData(generateDataReceivedMsg(this,workerData{i},i));
                    notify(this,"ActionMessageReceived",ed);
                    
                    % update the episode count
                    episodeCount = episodeCount + 1;
                    
                    % learn from worker data
                    learnFromWorkerData(this,workerData{i},activeSims);             
                    
                    % update the training manager
                    info.EpisodeInfo = workerEpisodeInfo{i};
                    info.SimulationInfo = simInfos{i};
                    info.WorkerID = i;
                    info.EpisodeCount = episodeCount;
                    
                    ed = rl.util.RLEventData(info);
                    notify(this,"EpisodeInfoReceivedFromWorker",ed);

                    % decrement the remaining sims
                    remainingSims = remainingSims - 1;
                end
                isFirstLoop = false;
            end
        end
    end
end