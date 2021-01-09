%% plot training graphs
close all; clear
agentN =  load('Agent_RT14N_13000.mat');

agentZ =  load('Agent_RT14Z_17100.mat');
% figure(1);clf;hold on;
% plot(movmean(agentZ.savedAgentResultStruct.TrainingStats.AverageReward,40));
% legend('Safety on',' Safety off')
figure('Name','Training Plot','Position',[1 1 600 300]) ; clf;hold on;
epsrewZ = movmean(agentZ.savedAgentResultStruct.TrainingStats.AverageReward,200);
epsrewN = movmean(agentN.savedAgentResultStruct.TrainingStats.AverageReward,200);
stop_eps =13000;
plot(epsrewZ(1:stop_eps),'LineWidth',2,'Color',[0 220/255 0]);
plot(epsrewN(1:stop_eps),'LineWidth',2,'Color',[1 51/255 51/255]);
xlim([0 stop_eps])
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
print(h,'trainingPlot.pdf','-dpdf','-r0')

%% tally up two experiences
num_iter = 500;

exp1=load("exp_Agent_N_6420.mat");
exp1= exp1.exp;
exp2=load("exp_Agent_Z_new_proprosed.mat");
exp2= exp2.exp;

[num1_0, num1_1, num1_2, num1_3, num1_4, num1_5, success1, avg1] = tally_up(exp1,num_iter);
[num2_0, num2_1, num2_2, num2_3, num2_4, num2_5, success2, avg2] = tally_up(exp2,num_iter);
both_success = success1 & success2;
rew1 = 0;
rew2 = 0;
for i = 1:500
if both_success(i) == 1
    rew1 = rew1 + sum( exp1(i).Reward.Data);
    rew2 = rew2 + sum( exp2(i).Reward.Data);
end
end
rew1/sum(both_success)
rew2/sum(both_success)
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