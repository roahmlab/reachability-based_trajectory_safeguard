%%

%%

w = warning ('off','all');
rng(0)
run_highway_testing

S.eval = 1; % if eval save episode as file
S.plot_sim_flag = 0;
AH.plot_flag = 0; % for paper plotting, usually keep off
S.plot_AH_flag = 0; % AH plotting, toggle for reference, turn off when AH.plot_flag is on;
S.epscur = 1;%61;this is for video generating episode

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


%%
S.safety_layer = 'R';
% load('Agent_Z_RT8_24300.mat');
%  load('Agent_Z_RT7_38166.mat');% 82%goal
%  load('Agent_N_RT8_28700.mat');
%  load('savedAgents_RT10ZD/Agent6800.mat')
  load('Agent_RT14N_12200.mat')

tic
simOptions = rlSimulationOptions('MaxSteps',200,'NumSimulations',500,'UseParallel',0);
exp = sim(env,saved_agent,simOptions);
toc
save("exp_Agent_RT14N_12200.mat","exp");
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run_highway_testing
S.eval = 1;
S.plot_sim_flag = 1;
%% Setup Environment
obsInfo = rlNumericSpec([5 1]); % 4*6 + 3 % 3*6+3
obsInfo.Name = 'Highway States';

%first arg spd, second lane y
actInfo = rlNumericSpec([2 1],'LowerLimit',[-1; -1],'UpperLimit',[1 ;1]);
actInfo.Name = 'act';

stepfun = @(action,LoggedSig)S.step(action);
reset_fun = @()S.reset();
env = rlFunctionEnv(obsInfo,actInfo,stepfun,reset_fun);


w = warning ('off','all');
run_highway_testing
S.eval = 1;
S.plot_sim_flag = 1;
AH.plot_flag = 0;
S.epscur =61;
%% Setup Environment
obsInfo = rlNumericSpec([5 1]); % 4*6 + 3 % 3*6+3
obsInfo.Name = 'Highway States';

%first arg spd, second lane y
actInfo = rlNumericSpec([2 1],'LowerLimit',[-1; -1],'UpperLimit',[1 ;1]);
actInfo.Name = 'act';

stepfun = @(action,LoggedSig)S.step(action);
reset_fun = @()S.reset();
env = rlFunctionEnv(obsInfo,actInfo,stepfun,reset_fun);


%%
S.safety_layer = 'N';
 load('Agent_N_RT8_35000.mat');

tic
simOptions = rlSimulationOptions('MaxSteps',200,'NumSimulations',1,'UseParallel',0);
exp = sim(env,saved_agent,simOptions);
toc
save("exp_Agent_N_RT9_35000.mat","exp");
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run_highway_testing
S.eval = 1;
S.plot_sim_flag = 0;
%% Setup Environment
obsInfo = rlNumericSpec([5 1]); % 4*6 + 3 % 3*6+3
obsInfo.Name = 'Highway States';

%first arg spd, second lane y
actInfo = rlNumericSpec([2 1],'LowerLimit',[-1; -1],'UpperLimit',[1 ;1]);
actInfo.Name = 'act';

stepfun = @(action,LoggedSig)S.step(action);
reset_fun = @()S.reset();
env = rlFunctionEnv(obsInfo,actInfo,stepfun,reset_fun);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
num_iter= 500;
num_5 = 0;
num_4 = 0;
num_3 = 0;
num_2 = 0;
num_1 = 0;
num_0 = 0;
  
for i = 1:num_iter
    data = exp(i).IsDone.Data(end);
    if sum(abs((data - 2))<0.1)
        a = 1;
    end
    
    num_5 = num_5 + sum(abs((data - 5))<0.1);
    num_4 = num_4 + sum(abs((data - 4))<0.1);
    num_3 = num_3 + sum(abs((data - 3))<0.1);
    num_2 = num_2 + sum(abs((data - 2))<0.1);
    num_1 = num_1 + sum(abs((data - 1))<0.1);
    num_0 = num_0 + sum(abs((data - 0))<0.1);
end
%%
num_iter= 500;
total_reward = 0;
episode_count = 0;
for i = 1:num_iter
    if exp(i).IsDone.Data(end) == 5
    num_steps = size(exp(i).Reward.Data,1);
    
    total_reward = total_reward + sum(exp(i).Reward.Data);
    episode_count= episode_count +1;
    end
   
end
%%
num_iter= 500;
res = zeros(num_iter,num_steps);
for i = 1:num_iter
    howmany = exp(i).Observation.HighwayStates.Data(end-1,:,:);
    for j = 1:length(howmany(:))
       res(i,j)= any([2 6 10]-exp(i).Observation.HighwayStates.Data(end-1,:,j) < 0.3);
    end
end

nnz(res)/num_iter/num_steps
%%
num_5
num_4
num_3
num_2
num_1
num_0
episode_count
total_reward/num_0
% %% simulate two runs
% % simulate_saved_file(LoggedSignals)
% LoggedSignals=load('sim_summary_N_5-157_17-04-51.504.mat');
% LoggedSignals=LoggedSignals.LoggedSignals;
% LoggedSignals2=load('sim_summary_N_5-eval_20.mat');
% LoggedSignals2= LoggedSignals2.LoggedSignals;
% 
% W = dynamic_car_world('bounds',bounds,'N_obstacles',N_obstacles,'buffer',world_buffer,...
%     'verbose',verbose_level,'goal_radius',goal_radius,...
%     'obstacle_size_bounds',obstacle_size_bounds) ;
% A = highway_cruising_agent;
% A.state = LoggedSignals.summary.agent_info.state;
% A.time =  LoggedSignals.summary.agent_info.time;
% W.setup(LoggedSignals.envCars_hist{1})
% 
% W2 = dynamic_car_world('bounds',bounds,'N_obstacles',N_obstacles,'buffer',world_buffer,...
%     'verbose',verbose_level,'goal_radius',goal_radius,...
%     'obstacle_size_bounds',obstacle_size_bounds) ;
% A2 = highway_cruising_agent;
% A2.state = LoggedSignals2.summary.agent_info.state;
% A2.time =  LoggedSignals2.summary.agent_info.time;
% W2.setup(LoggedSignals2.envCars_hist{1})
% 
% figure(1);clf;hold on; axis equal;
% 
% for i = 1: A.time(end)/2
%     if ~isempty(LoggedSignals.envCars_hist{i})
%         subplot(2,1,1);cla;hold on; axis equal;
%         W.plot()
%         A.plot_at_time(i*2)
%         text(LoggedSignals.envCars_hist{i}(1,1) + 30, 6,num2str(LoggedSignals.envCars_hist{i}(1,2)))
%         xlim([LoggedSignals.envCars_hist{i}(1,1)-S.car_max_vis_dist, LoggedSignals.envCars_hist{i}(1,1)+S.car_max_vis_dist])
%     end  
%     if ~isempty(LoggedSignals2.envCars_hist{i})
%         subplot(2,1,2);cla;hold on; axis equal;
%         W2.plot()
%         A2.plot_at_time(i*2)
%         text(LoggedSignals2.envCars_hist{i}(1,1) + 30, 6,num2str(LoggedSignals2.envCars_hist{i}(1,2)))
%         xlim([LoggedSignals2.envCars_hist{i}(1,1)-S.car_max_vis_dist, LoggedSignals2.envCars_hist{i}(1,1)+S.car_max_vis_dist]) 
%     end
%     drawnow
%     %  pause(0.2)
% end





%%
function [num_0, num_1, num_2,num_3, num_4, num_5,success_trials]=tally_up(exp,num_iter)
num_5 = 0;
num_4 = 0;
num_3 = 0;
num_2 = 0;
num_1 = 0;
num_0 = 0;
success_trials = zeros(num_iter,1);
for i = 1:num_iter
    data = exp(i).IsDone.Data;
    
    num_5 = num_5 + sum(abs((data - 5))<0.1);
    num_4 = num_4 + sum(abs((data - 4))<0.1);
    num_3 = num_3 + sum(abs((data - 3))<0.1);
    num_2 = num_2 + sum(abs((data - 2))<0.1);
    num_1 = num_1 + sum(abs((data - 1))<0.1);
    if sum(abs((data - 5))<0.1) ~= 0
        success_trials (i) = 1;
    end
    num_0 = num_0 + sum(abs((data - 0))<0.1);
end
end

