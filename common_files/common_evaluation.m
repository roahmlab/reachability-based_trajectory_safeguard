%% plot training graphs
close all; clear
N_color = [252,141,98]/255;
ZD_color = [100,130,255]/255;
Z_color = ([102,250,140]-10)/255; 
%% car
plot_N = true;
plot_ZD = true;
agentN =  load('Agent_RT14N_35000.mat');
agentZ =  load('Agent_RT15ZNoise_21737.mat');
agentZD = load('Agent_RT14ZD_35000.mat');
big_avg = 600;
window_size = 600;
stop_eps =20000;
small_avg = 30;
%% drone
plot_N = false;
plot_ZD = true;
agentN =  load('Agent_D7Super2_N2988.mat');% keep crashing in very few steps, not training
agentZD =  load('Agent_D7ZDSuper2_2414.mat');
agentZ =  load('Agent_D7ZSuper2_4625.mat');
big_avg = 50;
window_size = 50;
stop_eps =1500;

%% cartpole
plot_N = true;
plot_ZD = false;
agentN =  load('Agent_CN_300.mat');
agentZ =  load('Agent300_Z.mat');
big_avg = 30;
window_size = 30;
stop_eps =300;
%%
% figure(1);clf;hold on;
% plot(movmean(agentZ.savedAgentResultStruct.TrainingStats.AverageReward,40));
% legend('Safety on',' Safety off')
figure('Name','Training Plot','Position',[1 1 600 300]) ; clf;hold on;
plot_minmax(agentZ,big_avg, window_size,stop_eps,Z_color)
plot_minmax(agentN,big_avg, window_size,stop_eps,N_color)
plot_minmax(agentZD,big_avg, window_size,stop_eps,ZD_color)

xlim([0 stop_eps])
% ylim([-80 131])
% legend('Safety On','Safety Off');
% legend('Location','southeast');
xlabel('Training Episode')
ylabel('Episode Reward')
set(gca,'FontSize',15); grid on
h = gcf;
ax= gca;
ax.XRuler.Exponent = 0;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'trainingPlot.pdf','-painters' ,'-dpdf','-r0')
%%
%% tally up car experience and compare
num_iter = 500;
only_success_trials = false;

exp1=load("exp_Agent_RT15ZNoise_19240.mat");
exp1= exp1.exp;
exp2=load("exp_Agent_RT14ZD_17820.mat");
exp2= exp2.exp;
exp3=load("exp_Agent_RT14N_17890.mat");
exp3= exp3.exp;
exp4=load("exp_Agent_RT14R_nohlp.mat");
exp4= exp4.exp;
exp5=load("exp_Agent_RT14R.mat");
exp5= exp5.exp;
% 0: executed as is outputed by RL
% 1; executed and caused a collision but didn't go through safety check
% 2: RL output edited by safety layer and was safe
% 3: RL output edited by safety layer but still was unsafe % should be 0
% 4: No replacement action found and executing failsafe action/ all 0 command
% 5: task completed
[num1_0, num1_1, num1_2, num1_3, num1_4, num1_5, success1, avg1] = tally_up(exp1,num_iter);
[num2_0, num2_1, num2_2, num2_3, num2_4, num2_5, success2, avg2] = tally_up(exp2,num_iter);
[num3_0, num3_1, num3_2, num3_3, num3_4, num3_5, success3, avg3] = tally_up(exp3,num_iter);
[num4_0, num4_1, num4_2, num4_3, num4_4, num4_5, success4, avg4] = tally_up(exp4,num_iter);
[num5_0, num5_1, num5_2, num5_3, num5_4, num5_5, success5, avg5] = tally_up(exp5,num_iter);

% 5 here should not be counted in for success
all_success = success1 & success2 & success3 & success4 & success5;
rew1 = [];
rew2 = [];
rew3 = [];
rew4 = [];
rew5 = [];
for i = 1:500 
if all_success(i) == 1 || ~only_success_trials % when things are stuck, simulation ends early, for that to not happen: intropolate the reward gained from last step to all remaining steps when calculating reward
    if ~only_success_trials && exp1(i).IsDone.Data(end) == 4
        rew1 = [rew1 (200-length(exp1(i).IsDone.Data))*exp1(i).Reward.Data(end)+sum(exp1(i).Reward.Data)];
    else
        rew1 = [rew1  sum(exp1(i).Reward.Data)];
    end % 200 here is the number of evaluation steps allowed
   
    if ~only_success_trials && exp2(i).IsDone.Data(end) == 4
        rew2 = [rew2 (200-length(exp2(i).IsDone.Data))*exp2(i).Reward.Data(end)+sum(exp2(i).Reward.Data)];
    else
        rew2 = [rew2  sum(exp2(i).Reward.Data)];
    end % 200 here is the number of evaluation steps allowedif ~only_success_trials && exp1(i).IsDone.Data(end) == 4
    
    if ~only_success_trials && exp3(i).IsDone.Data(end) == 4
        rew3 = [rew3 (200-length(exp3(i).IsDone.Data))*exp3(i).Reward.Data(end)+sum(exp3(i).Reward.Data)];
    else
        rew3 = [rew3  sum(exp3(i).Reward.Data)];
    end % 200 here is the number of evaluation steps allowedif ~only_success_trials && exp1(i).IsDone.Data(end) == 4
    if ~only_success_trials && exp4(i).IsDone.Data(end) == 4 
        rew4 = [rew4 (200-length(exp4(i).IsDone.Data))*exp4(i).Reward.Data(end)+sum(exp4(i).Reward.Data)];
    else
        rew4 = [rew4  sum(exp4(i).Reward.Data)];
    end % 200 here is the number of evaluation steps allowedif ~only_success_trials && exp1(i).IsDone.Data(end) == 4
    if ~only_success_trials && exp5(i).IsDone.Data(end) == 4
        rew5 = [rew5 (200-length(exp5(i).IsDone.Data))*exp5(i).Reward.Data(end)+sum(exp5(i).Reward.Data)];
    else
        rew5 = [rew5 sum(exp5(i).Reward.Data)];
    end % 200 here is the number of evaluation steps allowed
end
end
mean(rew1)
std(rew1)
mean(rew2)
std(rew2)
mean(rew3)
std(rew3)
mean(rew4)
std(rew4)
mean(rew5)
std(rew5)

    
% rew5/sum(all_success)

%%
num_iter = 500;
only_success_trials = false;
exp1=load("exp_Agent_D7Super2_Z1026.mat");
exp1= exp1.exp;
exp2=load("exp_Agent_D7Super2_ZD857.mat");
exp2= exp2.exp;
exp3=load("exp_Agent_D7Super2_N1203.mat");
exp3= exp3.exp;
exp4=load("exp_Agent_D7Super2_RNaive.mat");
exp4= exp4.exp;
exp5=load("exp_Agent_D7Super2_R.mat");
exp5= exp5.exp;
[num1_0, num1_1, num1_2, num1_3, num1_4, num1_5, success1, avg1] = tally_up(exp1,num_iter);
[num2_0, num2_1, num2_2, num2_3, num2_4, num2_5, success2, avg2] = tally_up(exp2,num_iter);
[num3_0, num3_1, num3_2, num3_3, num3_4, num3_5, success3, avg3] = tally_up(exp3,num_iter);
[num4_0, num4_1, num4_2, num4_3, num4_4, num4_5, success4, avg4] = tally_up(exp4,num_iter);
[num5_0, num5_1, num5_2, num5_3, num5_4, num5_5, success5, avg5] = tally_up(exp5,num_iter);
all_success = success1  & success2& success3 & success4 ; % 2 have too few successes
rew1 = [];
rew2 = [];
rew3 = [];
rew4 = [];
rew5 = [];
for i = 1:500 
if all_success(i) == 1 || ~only_success_trials % when things are stuck, simulation ends early, for that to not happen: intropolate the reward gained from last step to all remaining steps when calculating reward
    if ~only_success_trials && exp1(i).IsDone.Data(end) == 4
        rew1 = [rew1 (80-length(exp1(i).IsDone.Data))*exp1(i).Reward.Data(end)+sum(exp1(i).Reward.Data)];
    else
        rew1 = [rew1  sum(exp1(i).Reward.Data)];
    end % 200 here is the number of evaluation steps allowed
   
    if ~only_success_trials && exp2(i).IsDone.Data(end) == 4
        rew2 = [rew2 (80-length(exp2(i).IsDone.Data))*exp2(i).Reward.Data(end)+sum(exp2(i).Reward.Data)];
    else
        rew2 = [rew2  sum(exp2(i).Reward.Data)];
    end % 200 here is the number of evaluation steps allowedif ~only_success_trials && exp1(i).IsDone.Data(end) == 4
    
    if ~only_success_trials && exp3(i).IsDone.Data(end) == 4
        rew3 = [rew3 (80-length(exp3(i).IsDone.Data))*exp3(i).Reward.Data(end)+sum(exp3(i).Reward.Data)];
    else
        rew3 = [rew3  sum(exp3(i).Reward.Data)];
    end % 200 here is the number of evaluation steps allowedif ~only_success_trials && exp1(i).IsDone.Data(end) == 4
    if ~only_success_trials && exp4(i).IsDone.Data(end) == 4 
        rew4 = [rew4 (80-length(exp4(i).IsDone.Data))*exp4(i).Reward.Data(end)+sum(exp4(i).Reward.Data)];
    else
        rew4 = [rew4  sum(exp4(i).Reward.Data)];
    end % 200 here is the number of evaluation steps allowedif ~only_success_trials && exp1(i).IsDone.Data(end) == 4
    if ~only_success_trials && exp5(i).IsDone.Data(end) == 4
        rew5 = [rew5 (80-length(exp5(i).IsDone.Data))*exp5(i).Reward.Data(end)+sum(exp5(i).Reward.Data)];
    else
        rew5 = [rew5 sum(exp5(i).Reward.Data)];
    end % 200 here is the number of evaluation steps allowed
end
end
mean(rew1)
std(rew1)
mean(rew2)
std(rew2)
mean(rew3)
std(rew3)
mean(rew4)
std(rew4)
mean(rew5)
std(rew5)

    
%% tally up cartpole
num_iter = 500;
only_success_trials = false;
exp1=load("exp_Agent_Cartpole_RNaive.mat");
exp1= exp1.exp;
exp2=load("exp_Agent_Zrms22.mat");
exp2= exp2.exp;
exp3=load("exp_Agent_Nrms.mat");
exp3= exp3.exp;
[num1_0, num1_1, num1_2, num1_3, num1_4, num1_5, success1, avg1] = tally_up(exp1,num_iter);
[num2_0, num2_1, num2_2, num2_3, num2_4, num2_5, success2, avg2] = tally_up(exp2,num_iter);
[num3_0, num3_1, num3_2, num3_3, num3_4, num3_5, success3, avg3] = tally_up(exp3,num_iter);

rew1 = [];
rew2 = [];
rew3 = [];
limit = 4;%deg
success1= [];
success2= [];
success3= [];
for i = 1:num_iter
    success1 = [success1 determine_cartpole_success(exp1(i),limit)];
    success2 = [success2 determine_cartpole_success(exp2(i),limit)];
    success3 = [success3 determine_cartpole_success(exp3(i),limit)];

    if ~only_success_trials && exp1(i).IsDone.Data(end) == 4
        rew1 = [rew1 (100-length(exp1(i).IsDone.Data))*exp1(i).Reward.Data(end)+sum(exp1(i).Reward.Data)];
    else
        rew1 = [rew1  sum(exp1(i).Reward.Data)];
    end % 200 here is the number of evaluation steps allowed

    if ~only_success_trials && exp2(i).IsDone.Data(end) == 4
        rew2 = [rew2 (100-length(exp2(i).IsDone.Data))*exp2(i).Reward.Data(end)+sum(exp2(i).Reward.Data)];
    else
        rew2 = [rew2  sum(exp2(i).Reward.Data)];
    end % 200 here is the number of evaluation steps allowedif ~only_success_trials && exp1(i).IsDone.Data(end) == 4

    if ~only_success_trials && exp3(i).IsDone.Data(end) == 4
        rew3 = [rew3 (100-length(exp3(i).IsDone.Data))*exp3(i).Reward.Data(end)+sum(exp3(i).Reward.Data)];
    else
        rew3 = [rew3  sum(exp3(i).Reward.Data)];
    end
end
mean(rew1)
std(rew1)
mean(rew2)
std(rew2)
mean(rew3)
std(rew3)
%%
function success = determine_cartpole_success(exp,limit)
sin_limit = sin(deg2rad(limit));
cos_limit = cos(deg2rad(limit));
last_states_sin = exp.Observation.obs1.Data(1,:,end-10:end);
last_states_cos = exp.Observation.obs1.Data(2,:,end-10:end);
if all(abs(last_states_sin(:))<sin_limit) && all(last_states_cos(:)>cos_limit)
    success = 1;
else
    success = 0;
end
end

%% function plot minmax
function plot_minmax(agentZ,movmeansize, window_size,stop_eps,color)
epsrewZ       = movmean(agentZ.savedAgentResultStruct.TrainingStats.AverageReward,movmeansize); 
%% This can also be used for plotting min/max shading, but doesn't show really well so omitted.
epsrewZstd = movstd(agentZ.savedAgentResultStruct.TrainingStats.AverageReward,window_size);
% epsrewZmin = movmin(agentZ.savedAgentResultStruct.TrainingStats.AverageReward,window_size);
x_data = 1:length(epsrewZ);
x2 = [x_data, fliplr(x_data)];
inBetween = [(epsrewZ-epsrewZstd)', fliplr((epsrewZ+epsrewZstd)')];
pf = fill(x2, inBetween, color);
pf.FaceAlpha = 0.2;
pf.EdgeAlpha = 0;
endZ  = min([stop_eps length(epsrewZ)]);
plot(epsrewZ(1:endZ),'LineWidth',2,'Color',color);
% pZminmax =plot(x_data(1:endZ*2),epsrewZminmax(1:endZ*2),'--','LineWidth',5,'Color',color);
% pZmin =plot(epsrewZmin(1:endZ),'--','LineWidth',1,'Color',Z_color);
end
%% function to tally up
function [num_0, num_1, num_2,num_3, num_4, num_5,success_trials,avg_rew]=tally_up(exp,num_iter)
num_5 = 0;
num_4 = 0;
num_3 = 0;
num_2 = 0;
num_1 = 0;
num_0 = 0;
reward_total = 0;
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
        reward_total= reward_total+sum(exp(i).Reward.Data);
    end
    num_0 = num_0 + sum(abs((data - 0))<0.1);
end
avg_rew = reward_total/sum(success_trials);
end
