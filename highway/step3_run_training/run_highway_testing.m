clear; close all;
% plot flag

plot_sim_flag = 1;
plot_AH_flag = 1;
AH_debug_flag = 0;
%% set up required objects for simulation
% agent
bounds = [-500,0,-0.7,12.7] ;
goal_radius = 10 ;
world_buffer = 1 ; % this is used to make start and goal locations that are not too close to the robot or boundary

% planner
buffer = 0.3 ; % [m] distance used to buffer obstacles for planning, must be larger than the agent's footprint
t_plan = 2.0 ; % if t_plan = t_move, then real time planning is enforced
t_move = 2.0 ;
t_failsafe_move = 4.0;
Ts = t_move;
Tf = Ts*150;

% verbose level for printing
verbose_level = 0;

% HLP choices:
% 1. Highway HLP: give waypoint same as reward function
% waypoint far enough on that lane, use fmincon to track the waypoint
HLP = highway_HLP('bounds',bounds,'verbose',verbose_level);

% % 2. RRT* HLP: use RRT* algorithm to proprose waypoints in the workspace
% % not inside the obstacles, use fmincon to track the waypoint
% RTD_HLP_grow_tree_mode ='seed';
% HLP_buffer = 0.1;
% HLP = RRT_star_HLP('timeout',0.5,'bounds',bounds,...                          %5m/s *2s peak time
%     'grow_tree_mode',RTD_HLP_grow_tree_mode,'new_node_max_distance_from_agent',5*2,'verbose',verbose_level);


% automated from here
A = highway_cruising_agent('plot_trajectory_at_time_flag', 0);
A.integrator_type= 'ode45';
A. desired_initial_condition = [bounds(1)+10; 6; 0; 22; 0];
W = dynamic_car_world('bounds',bounds,'buffer',world_buffer,...
    'verbose',verbose_level,'goal_radius',goal_radius) ;

AH = highwayAgentHelper(A,'zono_full_7.13_1spd.mat',HLP,'t_move',t_move,'t_failsafe_move',t_failsafe_move,'eps',0.001,'verbose',verbose_level,'plot_flag',AH_debug_flag);
S = rlsimulator(AH,W,'plot_sim_flag',plot_sim_flag, 'safety_layer','Z','plot_AH_flag',plot_AH_flag);
AH.S = S;
random_inputs = false; %true for user inputs using arrow keys, false for random inputs

S.eval = 0; %turn on evaluation so summary will be saved, and action will not have random components
% to the episode number

if true
    for j = 1:10
        S.reset();
        for i = 1: 150
            if random_inputs
                value = 114;
            else
                figure(1)
                k = waitforbuttonpress;
                value = double(get(gcf,'CurrentCharacter'));
            end
            
            % 28 leftarrow
            % 29 rightarrow
            % 30 uparrow
            % 31 downarrow
            
            if value == 28
                [~,~,IsDone,LoggedSignal]=S.step([-2/3;0]);
            elseif value == 29
                [~,~,IsDone,LoggedSignal]=S.step([1;0]);
            elseif value == 30
                [~,~,IsDone,LoggedSignal]=S.step([1/3;1]);
            elseif value == 31
                [~,~,IsDone,LoggedSignal]=S.step([1/3;-1]);
            elseif value == 114
                [~,~,IsDone,LoggedSignal]=S.step([rand*2-1;rand*2-1]);
            else
                [~,~,IsDone,LoggedSignal]=S.step([1/3;0]);
            end
            
            if IsDone == 1 || IsDone == 3 || IsDone == 4 || IsDone == 5
                %crash       
                %      crash with safety layer on
                %                      safely stopped but stuck
                %                                           reached goal!
                break
            end
            
        end
    end
end

done = 'Setup Complete'