function parCommTrain(trainMgr)
% PARCOMMTRAIN Uses the training manager to train an agent in parallel
% using the comm interface

% Revised: 10-26-2018
% Copyright 2018 The MathWorks, Inc.


% extract the agent, env, and training options from the train mgr
agent = trainMgr.Agent;
env = trainMgr.Environment;
trainingOptions = trainMgr.TrainingOptions;
maxEpisodes = trainingOptions.MaxEpisodes;
maxSteps = trainingOptions.MaxStepsPerEpisode;

pool = gcp();
numWorkers = pool.NumWorkers;

% use async or sync
useAsync = strcmpi(trainingOptions.Parallelization,"async");

% parse the parallel options
popts                = trainingOptions.ParallelizationOptions;
stepsUntilDataIsSent = popts.StepsUntilDataIsSent;
data2sendFromWorkers = popts.DataToSendFromWorkers;

% use blocking parameters? (workers wait until params are received)
useBlocking = popts.BlockWorkersUntilParametersAreReceived || ~useAsync;
% send gradient or experience to host
sendGradientsFlag = strcmpi(data2sendFromWorkers,'Gradients');

%% nested host functions
function nestedAsyncDataReceivedOnHostCallback(communicator,agent,message)
    % this is called when a worker sends experiences back to the host
    % ASYNC

    % get the worker id from the message
    wid = message.WorkerID;

    if ~agent.TerminateSimulation
        if sendGradientsFlag
            % get the gradients from the message
            gradientBuffer = message.Data;
            % update representations from gradients
            applyGradient(agent,gradientBuffer);
        else
            % get the experience from the message
            experiences = message.Data;
            % learn from the experience
            learnFromExperiences(agent,experiences);
        end
        % tell the world data has been received
        notifyDataReceivedFromWorker(trainMgr,message);
    end

    % get the parameters from the agent
    p = getLearnableParameters(agent);

    % send them back to the worker(s)
    if useBlocking
        sendParameters(communicator,p,wid);
    else
        % REVISIT this will cause MORE communication between the host and
        % workers which may slow down time/step. The benefit is that the
        % workers will have the latest parameter set SOONER than if params were
        % only send back to the worker that sent the data
        sendParameters(communicator,p);
    end

    drawnow();
end

receivedFromWorkersCount = 0;
function nestedSyncDataReceivedOnHostCallback(communicator,agent,message)
    % this is called when a worker sends experiences back to the host
    % SYNC

    if ~agent.TerminateSimulation
        % get the worker id from the message
        % wid = message.WorkerID;
        if sendGradientsFlag
            % get the gradients from the message
            gradientBuffer = message.Data;
            % update representations from gradients
            applyGradient(agent,gradientBuffer);
        else
            % get the experience from the message
            experiences = message.Data;
            % learn from the experience
            learnFromExperiences(agent,experiences);
        end
        % tell the world data has been received
        notifyDataReceivedFromWorker(trainMgr,message);
    end
    % record that the worker's exps were received
    receivedFromWorkersCount = receivedFromWorkersCount + 1;


    remainingEpisodes = maxEpisodes - trainMgr.EpisodeCount;
    activeSims = min(numWorkers,remainingEpisodes);
    if receivedFromWorkersCount >= activeSims
        % NOTE parameters must be sent back to the workers even if training was
        % cancelled since they will block on the worker until received

        % get the parameters from the agent
        p = getLearnableParameters(agent);

        % send the params back to all the workers
        sendParameters(communicator,p);

        % reset state
        receivedFromWorkersCount = 0;
    end
end

% REVISIT Ideally we would be able to place a listener on the environment
% for the "EpisodeFinished" event. As of now, fetchNext is too slow to
% get the simulation outputs when learning on the host.
finishedSimCount = 0;
startedSimCount = numWorkers;
continueTraining = true;
% siminfos = cell(maxEpisodes,1);
function nestedTrainingInfoReceivedCB(communicator,agent,message)
    % update the sim count
    finishedSimCount = finishedSimCount + 1;

    % get the worker id from the message
    wid = message.WorkerID;
    
    % get the sim info from the worker
    info = message.Data;
    info.WorkerID = wid;
    info.WorkerEpisodeCount = info.EpisodeCount;
    info.EpisodeCount = finishedSimCount;

%     siminfos{finishedSimCount} = message.Data.SimulationInfo;
    
    % did the policy request to stop?
    continueTraining = ~agent.TerminateSimulation;

    % update the training manager
    if continueTraining && finishedSimCount <= maxEpisodes
        stopTraining = update(trainMgr,info);

        % NOTE: If any of the following terms are FALSE, then
        % continueTraining will always be false, resulting in the training
        % being killed on all workers
        continueTraining = ...
            ~stopTraining && ...
            continueTraining;
    end

    % make sure we stop training if the max episodes has been
    % distributed (each worker is requested to sim N episodes)
    if continueTraining && startedSimCount < maxEpisodes

        % tell the workers to continue training
        sendContinueTraining(communicator,true,wid);
        startedSimCount = startedSimCount + 1;
    else
        % kill the simulation task
        sendContinueTraining(communicator,false,wid);
    end
end

%% create the host communicators
hostComm = rl.comm.ParHostCommunicator();

% update the experience callback depending upon mode
if useAsync
    datacb = @(communicator,agent,message) ...
        nestedAsyncDataReceivedOnHostCallback(communicator,agent,message);
else
    datacb = @(communicator,agent,message) ...
        nestedSyncDataReceivedOnHostCallback(communicator,agent,message);
end

% register the callbacks to the host communicator
if sendGradientsFlag
    registerGradientReceivedCallback(hostComm,agent,datacb);
else
    registerExperienceReceivedCallback(hostComm,agent,datacb);
end

infcb = @(communicator,agent,message) ...
    nestedTrainingInfoReceivedCB(communicator,agent,message);
registerTrainingInfoReceivedCallback(hostComm,agent,infcb);

%% create a commPolicy
% check for "send experiences at end of episode" mode
if stepsUntilDataIsSent == -1
    stepsUntilDataIsSent   = maxSteps;
    sendAtEndOfEpisodeMode = true;
else
    sendAtEndOfEpisodeMode = false;
end

% For PPO, always send data at the end of the episode
if isa(agent, 'rl.agent.rlPPOAgent')
    sendAtEndOfEpisodeMode = true;
end

% create the callbacks and make sure the agent is configured to use
% getActionWithExploration
if strcmpi(data2sendFromWorkers,'Gradients')
    setStepMode(agent,"sim-with-exploration-and-append");
    expfcn = @(agent,communicator,experience) localTransmitGradients(...
        agent,communicator,experience,stepsUntilDataIsSent);
else
    setStepMode(agent,"sim-with-exploration");
    expfcn = @(agent,communicator,experience) localTransmitExperience(...
        agent,communicator,experience,stepsUntilDataIsSent,sendAtEndOfEpisodeMode,useBlocking);
end

% create the policy object
preSimFcn  = @(agent,comm,simCount) localPreSimFcn   (agent,comm,simCount);
postSimFcn = @(agent,comm,info    ) localPostSimFcn  (agent,comm,info    );
setupFcn   = @(agent,comm         ) localSetupWorkers(agent,comm         );
commPolicy = rl.comm.CommPolicy(agent,hostComm,expfcn,preSimFcn,postSimFcn,setupFcn);

%% distribute the simulations
function nestedSimDistCB(comm)
    sendContinueTraining(comm,true);
    setActionMessage(trainMgr,getString(message('rl:general:TrainingManagerRunningTasksOnWorkers')));
end

% Add a listener to manually kill training (this likely isn't needed in the
% case where the host is doing the learning and parfeval is being called to
% manage the sims on the workers since they are cancelled directly)
trainlist(1) = event.listener(trainMgr,'TrainingManuallyStopped',...
    @(src,ed) sendManuallyStopTraining(hostComm));

% Add a listener to start sims once all the tasks have been distributed
trainlist(2) = event.listener(env,'SimulationsDistributed',...
    @(src,ed) nestedSimDistCB(hostComm)); 
clnuplist = onCleanup(@() delete(trainlist));

% set the action msg on the episode mgr
setActionMessage(trainMgr,getString(message('rl:general:TrainingManagerRunningTasksOnWorkers')));
clnupActMsg = onCleanup(@() setActionMessage(trainMgr,getString(message('rl:general:TrainingManagerCleaningUpWorkers'))));

simOptions = getSimulationOptions(trainingOptions);
simWithPolicy(env,commPolicy,simOptions);

% % attach the simulation info to the training manager
% trainMgr.SimulationInfo = vertcat(siminfos{:});
end

%% Local function executed on the workers
function localSetupWorkers(agent,comm)
% clear persistent variables in local functions that run on the workers
localTransmitExperience();
localTransmitGradients();

% register a training stopped cb
registerManuallyStoppedTrainingCallback(comm,agent,...
    @(comm,agent,msg) terminateSimulation(agent,true));
end

function localPreSimFcn(agent,communicator,simCount)
% don't start simulating until a connection is established with the host
if simCount < 2
    continueTraining = receiveContinueTraining(communicator);
    terminateSimulation(agent,~continueTraining);
end
end

function localPostSimFcn(agent,communicator,info) 
continueTraining = handshakeTrainingInfo(communicator,info);
terminateSimulation(agent,~continueTraining);
end

function localTransmitExperience(agent,communicator,experience,...
    stepsUntilExperiencesAreSent,sendAtEndOfEpisodeMode,useBlocking)
% called every time an experience is received on the worker (during the
% call to step)

persistent stepCount2SendExperiences;
persistent experienceBuffer;

if nargin
    if isempty(stepCount2SendExperiences)
        stepCount2SendExperiences = 0;
    end
    if isempty(experienceBuffer)
        experienceBuffer = cell(stepsUntilExperiencesAreSent,1);
    end
    
    % check to see if it is time to send experiences
    stepCount2SendExperiences = stepCount2SendExperiences + 1;
    experienceBuffer{stepCount2SendExperiences,1} = experience;
    if (stepCount2SendExperiences >= stepsUntilExperiencesAreSent) || ...
            (sendAtEndOfEpisodeMode && experience{5} > 0)
        
        % remove empty experiences
        experienceBuffer((stepCount2SendExperiences+1):end) = [];
        experienceBuffer = preprocessExperience(agent, experienceBuffer);
        
        if useBlocking
            % send the experience back to the host and get parameters back
            params = handshakeExperience(communicator,experienceBuffer);

            % install the parameters from the host
            setLearnableParameters(agent,params);
        else
            % send experiences back to host. parametes will be installed
            % via the registered localParametersReceiveOnWorker callback
            sendExperience(communicator,experienceBuffer);
        end

        % reset the step count and experience buffer
        stepCount2SendExperiences = 0;
        experienceBuffer = cell(stepsUntilExperiencesAreSent,1);
    end
else
    % clear persistent vars
    stepCount2SendExperiences = [];
    experienceBuffer = [];
end
end

function localTransmitGradients(agent,communicator,experience,...
    stepsUntilDataIsSent)
% called every time an experience is received on the worker (during the
% call to step)

persistent stepCount2SendData;
if nargin
    if isempty(stepCount2SendData)
        stepCount2SendData = 0;
    end
    
    % check to see if it is time to send experiences
    stepCount2SendData = stepCount2SendData + 1;
    % send every stepsUntilDataIsSent or at end of episode
    if (stepCount2SendData >= stepsUntilDataIsSent) || ...
            (experience{5} > 0)
        
        % accumulate gradients from agent's experience buffer
        gradientBuffer = accumulateGradient(agent,stepCount2SendData);
        
        % send the gradient to the host and get parameters back
        params = handshakeGradient(communicator,gradientBuffer);
        
        % install the parameters from the host
        setLearnableParameters(agent,params);

        % reset the step count and experience buffer
        stepCount2SendData = 0;
    end
else
    % clear persistent vars
    stepCount2SendData = [];
end
end
