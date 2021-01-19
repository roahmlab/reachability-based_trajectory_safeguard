run_cartpole_testing
%% Setup Environment
obsInfo = rlNumericSpec([5 1]);
numObservations = obsInfo.Dimension(1);
actInfo =rlNumericSpec(1) ;
actInfo.LowerLimit = -1;
actInfo.UpperLimit =  1;
reset_fun = @()Env.reset();
stepfun = @(action,LoggedSig)Env.step(action);


env = rlFunctionEnvNew(obsInfo,actInfo,stepfun,reset_fun);

% plotter = CartPoleVisualizerLimitedTrackNew(env,A);
% plotter.plot();

%% Create DDPG agent
statePath = [
    imageInputLayer([numObservations 1 1],'Normalization','none','Name','observation')
    fullyConnectedLayer(128,'Name','CriticStateFC1')
    reluLayer('Name','CriticRelu1')
    fullyConnectedLayer(200,'Name','CriticStateFC2')];

actionPath = [
    imageInputLayer([1 1 1],'Normalization','none','Name','action')
    fullyConnectedLayer(200,'Name','CriticActionFC1','BiasLearnRateFactor',0)];

commonPath = [
    additionLayer(2,'Name','add')
    reluLayer('Name','CriticCommonRelu')
    fullyConnectedLayer(1,'Name','CriticOutput')];

criticNetwork = layerGraph(statePath);
criticNetwork = addLayers(criticNetwork,actionPath);
criticNetwork = addLayers(criticNetwork,commonPath);
    
criticNetwork = connectLayers(criticNetwork,'CriticStateFC2','add/in1');
criticNetwork = connectLayers(criticNetwork,'CriticActionFC1','add/in2');

%% critic Option

criticOptions = rlRepresentationOptions('LearnRate',1e-02,'GradientThreshold',1);
critic = rlQValueRepresentation(criticNetwork,obsInfo,actInfo,...
    'Observation',{'observation'},'Action',{'action'},criticOptions);


%% Agent Option

actorNetwork = [
    imageInputLayer([numObservations 1 1],'Normalization','none','Name','observation')
    fullyConnectedLayer(128,'Name','ActorFC1')
    reluLayer('Name','ActorRelu1')
    fullyConnectedLayer(200,'Name','ActorFC2')
    reluLayer('Name','ActorRelu2')
    fullyConnectedLayer(1,'Name','ActorFC3')
    tanhLayer('Name','ActorTanh1')
    
    scalingLayer('Name','ActorScaling','Scale',max(actInfo.UpperLimit))];

actorOptions = rlRepresentationOptions('LearnRate',5e-04,'GradientThreshold',1);

actor = rlDeterministicActorRepresentation(actorNetwork,obsInfo,actInfo,...
    'Observation',{'observation'},'Action',{'ActorScaling'},actorOptions);
    
agentOptions = rlDDPGAgentOptions(...
    'SampleTime',0.1,...		    'SampleTime',Ts,...
    'TargetSmoothFactor',1e-3,...		    'TargetSmoothFactor',1e-3,...
    'ExperienceBufferLength',1e6,...		    'ExperienceBufferLength',1e6,...
    'MiniBatchSize',128);


agentOptions.NoiseOptions.Variance = 0.2;
agentOptions.NoiseOptions.VarianceDecayRate = 2e-4;

%% Buffer
% agent = rlTD3Agent(actor,critic,agentOptions);
agent = rlDDPGAgent(actor,critic,agentOptions);
%  load savedAgents_cartpole_Z_rmsprop2/Agent435.mat
%% Train Agent
maxepisodes = 1000;
%maxsteps = ceil(Tf/Ts);
maxsteps = 200;

trainOpts = rlTrainingOptions(...
    'MaxEpisodes', maxepisodes, ...		    'MaxEpisodes',maxepisodes,...
    'MaxStepsPerEpisode', maxsteps, ...		    'MaxStepsPerEpisode',maxsteps,...
    'ScoreAveragingWindowLength',5,...		    'ScoreAveragingWindowLength',5,...
    'Verbose', true, ...		    'Verbose',false,...
    'Plots','training-progress',...		    'Plots','training-progress',...
    'StopTrainingCriteria','AverageReward',...		    'StopTrainingCriteria','AverageReward',...
    'StopTrainingValue',1950, ... %480		    'StopTrainingValue',-400,...
    'SaveAgentCriteria','EpisodeReward',...
    'SaveAgentValue',-100);

%plot(env);

trainOpts.UseParallel = false; % true
% trainOpts.ParallelizationOptions.Mode = "async";
% trainOpts.ParallelizationOptions.DataToSendFromWorkers = "Experiences";
% trainOpts.ParallelizationOptions.StepsUntilDataIsSent = 32;
trainingStats = train(agent,env,trainOpts);
%%
% load cartpole_agent_2000.mat
% %% Result
simOptions = rlSimulationOptions('MaxSteps',1000);


experience = sim(env,saved_agent,simOptions);
% totalReward = sum(experience.Reward);