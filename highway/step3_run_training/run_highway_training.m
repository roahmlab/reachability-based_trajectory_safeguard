run_highway_testing
S.plot_sim_flag = 1;
S.eval = 0;
S.discrete_flag = false;
%% Setup Environment
obsInfo = rlNumericSpec([6 1]);
obsInfo.Name = 'Highway States';

%first arg spd, second lane y
actInfo = rlNumericSpec([2 1],'LowerLimit',[-1; -1],'UpperLimit',[1 ;1]);
actInfo.Name = 'act';

stepfun = @(action,LoggedSig)S.step(action);
reset_fun = @()S.reset();
env = rlFunctionEnv(obsInfo,actInfo,stepfun,reset_fun);
%% Setup Actor Network
actorNetwork = [
    imageInputLayer([6 1],'Normalization','none','Name','observation')
    fullyConnectedLayer(30,'Name','ActorFC1','WeightsInitializer','he')
    reluLayer('Name','ActorRelu1')
    fullyConnectedLayer(60,'Name','ActorFC2','WeightsInitializer','he')
    reluLayer('Name','ActorRelu2')
    fullyConnectedLayer(60,'Name','ActorFC3','WeightsInitializer','he')
    reluLayer('Name','ActorRelu3')
    fullyConnectedLayer(2,'Name','ActorFC4')% last layer tanh so keep using xavier
    tanhLayer('Name','ActorTanh1')];
lgraph = layerGraph(actorNetwork);
actorOptions = rlRepresentationOptions('LearnRate',5e-04,'GradientThreshold',1);

actor = rlDeterministicActorRepresentation(lgraph,obsInfo,actInfo,...
    'Observation',{'observation'},'Action',{'ActorTanh1'},actorOptions);


%% Create Critic Network
statePath = [
    imageInputLayer([6 1],'Normalization','none','Name','observation')
    fullyConnectedLayer(30,'Name','CriticStateFC1','WeightsInitializer','he')
    reluLayer('Name','CriticRelu1')
    fullyConnectedLayer(60,'Name','CriticStateFC2','WeightsInitializer','he')];

actionPath = [
    imageInputLayer([2 1],'Normalization','none','Name','action')
    fullyConnectedLayer(60,'Name','CriticActionFC1','BiasLearnRateFactor',0,'WeightsInitializer','he')];

commonPath = [
    additionLayer(2,'Name','add')
    reluLayer('Name','CriticCommonRelu')
    fullyConnectedLayer(1,'Name','CriticOutput','WeightsInitializer','he')];

criticNetwork = layerGraph(statePath);
criticNetwork = addLayers(criticNetwork,actionPath);
criticNetwork = addLayers(criticNetwork,commonPath);

criticNetwork = connectLayers(criticNetwork,'CriticStateFC2','add/in1');
criticNetwork = connectLayers(criticNetwork,'CriticActionFC1','add/in2');

criticOptions = rlRepresentationOptions('LearnRate',1e-3,'GradientThreshold',1);
critic = rlQValueRepresentation(criticNetwork,obsInfo,actInfo,...
    'Observation',{'observation'},'Action',{'action'},criticOptions);


%% Agent Option
agentOptions = rlTD3AgentOptions(...
    'SampleTime',Ts,...
    'TargetSmoothFactor',1e-3,...
    'ExperienceBufferLength',1e6,...
    'MiniBatchSize',128);
agentOptions.ExplorationModel.Variance = [4; 4];
agentOptions.ExplorationModel.VarianceDecayRate = 1-0.9999992;
agentOptions.ExplorationModel.LowerLimit = [-2; -2];
agentOptions.ExplorationModel.UpperLimit = [2; 2];
agentOptions.ExplorationModel.Mean = [-2;-2];
rlagent = rlTD3Agent(actor,critic,agentOptions);

%% Training Agent
maxepisodes = 20000;
maxsteps = ceil(Tf/Ts);
trainOpts = rlTrainingOptions(...
    'MaxEpisodes', maxepisodes, ...
    'MaxStepsPerEpisode', 600, ...
    'Verbose', true, ...
    'Plots','training-progress',...
    'SaveAgentCriteria','EpisodeCount',...
    'SaveAgentValue',1,...
    'StopTrainingCriteria','AverageReward',...
    'StopTrainingValue',1e10);


trainOpts.UseParallel = ~S.plot_sim_flag;
trainOpts.ParallelizationOptions.Mode = "async";
trainOpts.ParallelizationOptions.DataToSendFromWorkers = "Experiences";
trainOpts.ParallelizationOptions.StepsUntilDataIsSent = 128;

trainingStats = train(rlagent,env,trainOpts);
