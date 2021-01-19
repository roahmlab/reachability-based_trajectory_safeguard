% clear;clc;close all;
% Author: Simon Shao
% Created: 7 Aug 2020
% Updated: 9 Aug 2020
% clear;clc;close all;
run_quadrotor_testing
S.plot_sim_flag = 0;
w = warning ('off','all');
%% Setup Environment
S.safety_layer = 'Z';
S.discrete_flag = 0;
Tf = 40;%s
Ts = S.AH.t_move;
num_ob  = 35;
num_act = 3;
obsInfo = rlNumericSpec([num_ob 1]); % 4*6 + 3 % 3*6+3
obsInfo.Name = 'Highway States';

%first arg spd, second lane y
actInfo = rlNumericSpec([num_act 1],'LowerLimit',-1*ones(num_act,1),'UpperLimit',1*ones(num_act,1));
actInfo.Name = 'act';

stepfun = @(action,LoggedSig)S.step(action);
reset_fun = @()S.reset();
env = rlFunctionEnv(obsInfo,actInfo,stepfun,reset_fun);
%% Setup Actor Network
actorNetwork = [
    imageInputLayer([num_ob 1],'Normalization','none','Name','observation')
    fullyConnectedLayer(100,'Name','ActorFC1','WeightsInitializer','he')
    reluLayer('Name','ActorRelu1')
    fullyConnectedLayer(200,'Name','ActorFC2','WeightsInitializer','he')
    reluLayer('Name','ActorRelu2')
    fullyConnectedLayer(200,'Name','ActorFC3','WeightsInitializer','he')
    reluLayer('Name','ActorRelu3')
    fullyConnectedLayer(num_act,'Name','ActorFC4')
    tanhLayer('Name','ActorTanh1')];
lgraph = layerGraph(actorNetwork);
actorOptions = rlRepresentationOptions('LearnRate',5e-04,'GradientThreshold',1);

actor = rlDeterministicActorRepresentation(lgraph,obsInfo,actInfo,...
    'Observation',{'observation'},'Action',{'ActorTanh1'},actorOptions);


%% Create Critic Network
statePath = [
    imageInputLayer([num_ob 1],'Normalization','none','Name','observation')
    fullyConnectedLayer(200,'Name','CriticStateFC1','WeightsInitializer','he')
    reluLayer('Name','CriticRelu1')
    fullyConnectedLayer(100,'Name','CriticStateFC2','WeightsInitializer','he')];

actionPath = [
    imageInputLayer([num_act 1],'Normalization','none','Name','action')
    fullyConnectedLayer(100,'Name','CriticActionFC1','BiasLearnRateFactor',0,'WeightsInitializer','he')];

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
agentOptions.ExplorationModel.Variance = 2*ones(num_act,1);
agentOptions.ExplorationModel.VarianceDecayRate = 1-0.9999992;
agentOptions.ExplorationModel.LowerLimit = -1*ones(num_act,1);
agentOptions.ExplorationModel.UpperLimit =  1*ones(num_act,1);
agentOptions.ExplorationModel.Mean = -1*ones(num_act,1);
rlagent = rlTD3Agent(actor,critic,agentOptions);
% load savedAgentsD7Z/Agent101.mat
%% Training Agent
maxepisodes = 40000;
maxsteps = ceil(Tf/Ts);
trainOpts = rlTrainingOptions(...
    'MaxEpisodes', maxepisodes, ...
    'MaxStepsPerEpisode', maxsteps, ...
    'Verbose', true, ...
    'Plots','training-progress',...
    'SaveAgentCriteria','EpisodeCount',...
    'SaveAgentValue',1,...
    'StopTrainingCriteria','AverageReward',...
    'StopTrainingValue',1e10,'StopOnError','off');


trainOpts.UseParallel = ~S.plot_sim_flag;
trainOpts.ParallelizationOptions.Mode = "async";
trainOpts.ParallelizationOptions.DataToSendFromWorkers = "Experiences";
trainOpts.ParallelizationOptions.StepsUntilDataIsSent = 128;

trainingStats = train(rlagent,env,trainOpts);
