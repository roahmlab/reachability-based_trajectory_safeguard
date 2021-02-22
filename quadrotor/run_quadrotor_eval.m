run_quadrotor_testing
S.eval = 1;
S.plot_AH_flag = 1;
S.plot_sim_flag = 1;
S.safety_layer = 'Z';
S.epscur =1;% plotting episode =  38;
S.discrete_flag = 0;
%%
w = warning ('off','all');

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
%%
% load('Agent_D7Super2_N1203.mat');
load('Agent_D7Super2_Z1026.mat');

tic
simOptions = rlSimulationOptions('MaxSteps',80,'NumSimulations',500,'UseParallel',0);
exp = sim(env,saved_agent,simOptions);
toc
save("exp_Agent_D7Super2_N1203.mat","exp");
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
