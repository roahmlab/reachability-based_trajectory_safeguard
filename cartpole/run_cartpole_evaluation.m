run_cartpole_testing
%% Setup Environment
AH.flags.safety_layer = 'R';

obsInfo = rlNumericSpec([5 1]);
% numObservations = obsInfo.Dimension(1);
actInfo =rlNumericSpec(1) ;
actInfo.LowerLimit = -1;
actInfo.UpperLimit =  1;
Env.evalflag = 1;
reset_fun = @()Env.reset();
stepfun = @(action,LoggedSig)Env.step(action);

env = rlFunctionEnv(obsInfo,actInfo,stepfun,reset_fun);

%
load('Agent163_Z_rms22.mat');
% load('Agent_CN_300.mat');
tic
simOptions = rlSimulationOptions('MaxSteps',100,'NumSimulations',500,'UseParallel',0);
exp = sim(env,saved_agent,simOptions);
toc
save("exp_Agent_Nrms.mat","exp");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run_cartpole_testing
AH.flags.safety_layer = 'Z';
%% Setup Environment
obsInfo = rlNumericSpec([5 1]);
numObservations = obsInfo.Dimension(1);
actInfo =rlNumericSpec(1) ;
actInfo.LowerLimit = -1;
actInfo.UpperLimit =  1;
Env.evalflag = 1;
reset_fun = @()Env.reset();
stepfun = @(action,LoggedSig)Env.step(action);


env = rlFunctionEnvNew(obsInfo,actInfo,stepfun,reset_fun);

%%
load('Agent300_N.mat');

tic
simOptions = rlSimulationOptions('MaxSteps',200,'NumSimulations',500,'UseParallel',0);
exp = sim(env,saved_agent,simOptions);
toc
% save("exp_Agent_N_3000.mat","exp");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
