classdef CommPolicy < rl.agent.AbstractPolicy
% COMMPOLICY

% Revised: 10-16-2018
% Copyright 2018 The MathWorks, Inc.

    properties
        
        % NOTE: All of the following properties must be serializable to the
        % workers. The exception is the WorkerCommunicator and listeners
        % which will be created when loaded on the worker
        
        % agent that will be wrapped
        Agent
        % host communicator
        HostCommunicator
        % transmitters
        PreSimTransmitter
        StepTransmitter
        PostSimTransmitter
        % setup comm on worker function
        SetupOnWorkerFcn
    end
    properties (Transient)
        % worker communicator
        WorkerCommunicator % built when loaded on the worker
        % sim listeners
        SimListeners
    end
    methods
        function this = CommPolicy(agent,hostcomm,stepTransmitter,preSimTransmitter,postSimTransmitter,setupFcn)
            this.Agent = agent;
            this.WorkerCommunicator = [];
            this.HostCommunicator = hostcomm;
            this.StepTransmitter = stepTransmitter;
            this.PreSimTransmitter = preSimTransmitter;
            this.PostSimTransmitter = postSimTransmitter;
            if nargin < 6
                setupFcn = [];
            end
            this.SetupOnWorkerFcn = setupFcn;
            buildListeners(this);
        end
        function delete(this)
            delete(this.SimListeners);
        end
        function ts = getSampleTime(this)
            ts = getSampleTime(this.Agent);
        end
        function info = getObservationInfo(this)
            info = getObservationInfo(this.Agent);
        end        
        function info = getActionInfo(this)
            info = getActionInfo(this.Agent);
        end
        function transmit(this,transmitter,data)
            if ~isempty(transmitter)
                try
                    transmitter(this.Agent,this.WorkerCommunicator,data);
                catch ex
                    rethrow(ex);
                end
            end
        end
    end
    methods (Hidden)
        function preSimFcn(this,simCount)
            
            % overload the postSimFcn to do a postSimTransmit
            preSimFcn@rl.agent.AbstractPolicy(this,simCount);
            
            % initialize the episode
            preSimFcn(this.Agent,simCount);
            % transmit the simcout at the start of sim
            transmit(this,this.PreSimTransmitter,simCount);
        end
        function postSimFcn(this,simCount,simInfo)
            % overload the postSimFcn to do a postSimTransmit
            postSimFcn@rl.agent.AbstractPolicy(this,simCount,simInfo);
            % cleanup the episode
            postSimFcn(this.Agent,simCount,simInfo);
            % transmit at the end of sim
            info.EpisodeInfo    = getEpisodeInfo(this);
            info.EpisodeCount   = simCount;
            info.SimulationInfo = simInfo;
            transmit(this,this.PostSimTransmitter,info);
        end
    end
    methods (Access = protected)
        function buildWorkerComm(this)
            % check if the worker comm is already built 
            if isempty(this.WorkerCommunicator)
                % build from the host comm
                this.WorkerCommunicator = buildParWorkerCommunicator(this.HostCommunicator);
                % register any other function via the SetupOnWorkerFcn
                if ~isempty(this.SetupOnWorkerFcn)
                    this.SetupOnWorkerFcn(this.Agent,this.WorkerCommunicator);
                end
            end
        end
        function buildListeners(this)
            this.SimListeners = addlistener(this.Agent,'TerminateSimulation','PostSet',...
                @(src,ed) terminateSimulation(this,this.Agent.TerminateSimulation));
        end
        function action = stepImpl(this,experience)
            
            % call step on the agent
            action = stepImpl(this.Agent,experience);
            
            % check for early termination REVISIT avoid duplicate call (by stepImpl on agent)
            experience = checkForEarlyTermination(this.Agent,experience);
            
            % transmit the experience REVISIT should this be before/after
            % the agent step? or maybe we should have a pre/post step
            % transmitter?
            transmit(this,this.StepTransmitter,experience);
        end
        function action = getActionImpl(this,observation)
            action = getActionImpl(this.Agent,observation);
        end
        function action = getInitialActionImpl(this,observation)
            action = getInitialActionImpl(this.Agent,observation);
        end
        function resetImpl(this)
            resetImpl(this.Agent);
        end
    end
    methods (Static)
        function this = loadobj(s)
            this = s;
            % build the worker comm
            buildWorkerComm(this);
            % build listeners
            buildListeners(this);
        end
    end
end