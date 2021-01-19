clear all;
A = dummy_cartpole_agent('integrator_type','ode45');
AH = cartpole_agentHelper(A,'cartpole_frs_9_22.mat','safety_layer','Z');
Env= cartpole_env(AH);
num_range = 5;
N_theta = num_range;% initial
N_theta_dot = num_range;% initial
N_x_dot = num_range;% initial & parameter
N_ka    = num_range;% initi2al & parameter
N_kpeak = num_range;% selectable &  parameter

theta_start = -pi; theta_end = pi;
theta_dot_start= -2; theta_dot_end = 2;
x_dot_start = -2; x_dot_end = 2;
ka_start = -1; ka_end = 1;
kpeak_start = -2; kpeak_end = 2;


% sample the parameter space
theta = linspace(theta_start,theta_end,N_theta);
theta_dot = linspace(theta_dot_start,theta_dot_end,N_theta_dot);
x_dot = linspace(x_dot_start,x_dot_end,N_x_dot);
ka = linspace(ka_start,ka_end,N_ka);
kpeak = linspace(kpeak_start,kpeak_end,N_kpeak);

% kv   = 0;
% ka   = 0;
y = [0 0 pi 0];

states = [];
t_all = [];
refU = [];
refZ = [];
refT = [];
t_now = 0;
max_v = 1;
% env = rlPredefinedEnvNew('CartPole-ContinuousLimitedTrack');
% simon = 0;

if false
%     plotter = rl.env.viz.CartPoleVisualizerLimitedTrack(env);
Env.reset();    
for i = 1:600
        k = waitforbuttonpress;
        value = double(get(gcf,'CurrentCharacter'));
        if value == 49
            k_pk  = -max_v;
        elseif value == 50
            k_pk =  -max_v/2;
        elseif value == 51
            k_pk = max_v/2;
        elseif value == 52
            k_pk = max_v;
        else
            k_pk = 0;
        end
        %     k_pk= rand()*max_v*2 -max_v;
         Env.step(k_pk);
%         if simon
%             regen = 1;
%             while regen
%                 [T,U,Z] =parameterized_cartpole_traj(k_pk,kv,ka);
%                 %     Z = [0 Z];
%                 %     T = [0 T];
%                 %     U = [U(1) U];
%                 Z= Z + y(end,1);
%             end
%          222   [t, y] = ode45(@(t,y)A.cartpole_dynamic_new(t,y,T,U,Z),[0 Ts],y(end,:));
%             kv = y(end,2); ka = (y(end,2)- y(end-1,2))/(t(end)-t(end-1));
%             states = [states;y];
%             refU = [refU U(1:100)];
%             refZ = [refZ Z(1:100)];
%             refT = [refT t_now+T(1:100)];
%             t_all =[t_all; t_now+t];
%             t_now = t_now + Ts;
%             for i = 1:size(t)
%                 env.State = y(i,:);
%                 plotter.plot();
%                 drawnow;3
%             end
%         end
    end
end

done = 'Setup Complete'