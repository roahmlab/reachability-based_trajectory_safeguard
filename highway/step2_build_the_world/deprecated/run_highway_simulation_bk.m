%% description
% This script runs a simulation with the TurtleBot in the simulator
% framework, using RRT to plan online.
%
% Author: Shreyas Kousik
% Created: 19 Oct 2019
% Updated: 24 Jan 2019
% Updated: 23 Feb 2020 Simon Shao 
% backup old
clear all;
addpath(genpath('../simulator/'));
addpath(genpath('../RTD/'));
addpath(genpath('../RTD_tutorial/'));
%%
plot_while_running = 1; 

%% user parameters
% agent
w = warning ('on','all');
desired_speed = 0.75 ; % m/s
sensor_radius = 2.5 ; % make this larger if initialize_tree_mode is 'once'

% world
obstacle_size_bounds = [1, 1] ; % side length [min, max]
N_obstacles = 0 ;
bounds = [-500,500,-0.7,12.7] ;
goal_radius = 10 ;
world_buffer = 1 ; % this is used to make start and goal locations that are not too close to the robot or boundary

% planner
buffer = 0.3 ; % [m] distance used to buffer obstacles for planning, must be larger than the agent's footprint
t_plan = 1.0 ; % if t_plan = t_move, then real time planning is enforced
t_move = 0.5 ;
Ts = t_move;
Tf = Ts*600;
initialize_tree_mode = 'iter' ; % 'iter' or 'once'
HLP_grow_tree_mode = 'new' ; % 'new' or 'seed' or 'keep' (only matters if using 'iter' above)
grow_tree_once_timeout = 10 ;
HLP_type = 'rrt' ; % 'rrt' or 'rrt*' or 'connect' or 'connect*'
new_node_growth_distance = 0.5 ;

% plotting
plot_HLP_flag = true;
plot_tree_growth_flag = true ;


% simulation
verbose_level = 0;
max_sim_iterations = 100 ;
max_sim_time = 100 ;

% automated from here
A = highway_cruising_agent;

% this is needed for the agent to track the RRT* output
A.LLC.gains.yaw = 10 ;
A.LLC.lookahead_time = 0.1 ;

W = dynamic_car_world('bounds',bounds,'N_obstacles',N_obstacles,'buffer',world_buffer,...
                     'verbose',verbose_level,'goal_radius',goal_radius,...
                     'obstacle_size_bounds',obstacle_size_bounds) ;
                 
P = zonotope_planner('t_move',t_move,'t_plan',t_plan);


%%
S = highway_sim_simulator(A,W,P,'allow_replan_errors',true,'verbose',verbose_level,...
              'max_sim_time',max_sim_time,'max_sim_iterations',max_sim_iterations,'plot_while_running',plot_while_running);
max_del = [];
max_h = [];
%%
 for j = 1:100
S.reset(1,1);
for i = 1: 200
     k = waitforbuttonpress;
     value = double(get(gcf,'CurrentCharacter'));
% if i == 10
%     value = 28;
% else
%  value = 5;
% end
% value =114;
    % 28 leftarrow
% 29 rightarrow
% 30 uparrow
% 31 downarrow
if value == 28
%     
    [~,~,IsDone,LoggedSignal]=S.step([-1;0],1,1);
elseif value == 29
    [~,~,IsDone,LoggedSignal]=S.step([1;0],1,1);
elseif value == 30
     [~,~,IsDone,LoggedSignal]=S.step([1/3;1],1,1);
 elseif value == 31
    [~,~,IsDone,LoggedSignal]=S.step([1/3;-1],1,1);
elseif value == 114
    [~,~,IsDone,LoggedSignal]=S.step([rand*2-1;rand*2-1],1,1);
else
     [~,~,IsDone,LoggedSignal]=S.step([1/3;0],1,1);
end



    if IsDone == 1 || IsDone == 3
        break
    end
%     max_del = [max_del; abs(S.A.state(5,end))];
%     max_h = [max_h; abs(S.A.state(3,end))];
%     if abs(S.A.state(5,end))>max_del
%         max_del = abs(S.A.state(5,end));
%     end
%     if abs(S.A.state(3,end))>max_h
%         max_h = abs(S.A.state(3,end));
%     end
end
 end
 %
 %%
%  figure(1)
% histogram(max_del);
% figure(2)
% histogram(max_h)
 %%

%% 
obsInfo = rlNumericSpec([21 1]); % 4*6 + 3 % 3*6+3
obsInfo.Name = 'Highway States';

%first arg spd, second lane y
actInfo = rlNumericSpec([2 1],'LowerLimit',[-1; -1],'UpperLimit',[1 ;1]);
actInfo.Name = 'act';

stepfun = @(action,LoggedSig)S.step(action,1,1);
reset_fun = @()S.reset(1,1);
env = rlFunctionEnv(obsInfo,actInfo,stepfun,reset_fun);
%%
actorNetwork = [
    imageInputLayer([21 1],'Normalization','none','Name','observation')
    fullyConnectedLayer(128,'Name','ActorFC1')
    reluLayer('Name','ActorRelu1')
    fullyConnectedLayer(200,'Name','ActorFC2')
    reluLayer('Name','ActorRelu2')
    fullyConnectedLayer(2,'Name','ActorFC3')
    tanhLayer('Name','ActorTanh1')];
lgraph = layerGraph(actorNetwork);
plot(lgraph)
actorOptions = rlRepresentationOptions('LearnRate',5e-04,'GradientThreshold',1);

actor = rlDeterministicActorRepresentation(lgraph,obsInfo,actInfo,...
    'Observation',{'observation'},'Action',{'ActorTanh1'},actorOptions);


%% Create DDPG agent
% critic_act_info = rlNumericSpec([1 1]);
statePath = [
    imageInputLayer([21 1],'Normalization','none','Name','observation')
    fullyConnectedLayer(128,'Name','CriticStateFC1')
    reluLayer('Name','CriticRelu1')
    fullyConnectedLayer(200,'Name','CriticStateFC2')];

actionPath = [
    imageInputLayer([2 1],'Normalization','none','Name','action')
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
plot(criticNetwork)
%% critic Option

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
agentOptions.ExplorationModel.VarianceDecayRate = 1-0.999997;
agentOptions.ExplorationModel.LowerLimit = [-2; -2];
agentOptions.ExplorationModel.UpperLimit = [2; 2];
agentOptions.ExplorationModel.Mean = [-2;-2];
agent = rlTD3Agent(actor,critic,agentOptions);
% % agentOptions.NoiseOptions.Variance = 0.4;

%% Train Agent
maxepisodes = 2000;
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


trainOpts.UseParallel = ~plot_while_running;
trainOpts.ParallelizationOptions.Mode = "async";
trainOpts.ParallelizationOptions.DataToSendFromWorkers = "Experiences";
trainOpts.ParallelizationOptions.StepsUntilDataIsSent = 128;
%%
trainingStats = train(agent,env,trainOpts);

%% Result
simOptions = rlSimulationOptions('MaxSteps',600,'NumSimulations',100,'UseParallel',0);
exp = sim(env,agent,simOptions);
% totalReward = sum(experience.Reward);
%% 
 %trainingStats = train(agent,env,trainingOpts);
 %%
  save("agent_after_static_obs_6.2.mat","agent","trainingStats")
 %%
  save("exp_after_more_lane_small_act_TD3_5.25_strategy2.mat","exp")
