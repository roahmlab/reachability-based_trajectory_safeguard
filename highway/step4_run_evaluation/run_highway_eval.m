w = warning ('off','all');
rng(0)
run_highway_testing

S.eval = 1; % if eval save episode as file
S.plot_sim_flag = 1;
AH.plot_flag = 1; % for paper plotting, usually keep off
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
S.safety_layer = 'Z';
load('Agent_RT15ZNoise_19240.mat');
%  load('Agent_D7Super2_ZD857.mat');
%  load('Agent_D7Super2_N1203.mat');
total_sim_tic = tic;
simOptions = rlSimulationOptions('MaxSteps',200,'NumSimulations',500,'UseParallel',0);
exp = sim(env,saved_agent,simOptions);
total_time = toc(total_sim_tic);
% save("exp_RT14ZD_8215.mat","exp");
