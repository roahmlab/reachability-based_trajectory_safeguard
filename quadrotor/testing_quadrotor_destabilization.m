%% description
% This script tests the quadrotor_agent and quadrotor_agentHelper to
% attempt to destabilize the drone with its default Mellinger LLC
% controller. You can destabilize the drone pretty easily by setting the
% a_max variable to 5 or so.
%
% Author: Shreyas Kousik
% Created: 10 Aug 2020
% Updated: -
%
%% user parameters
% quadrotor parameters
v_max = 5 ; % m/s (don't change this)
a_max = 3 ; % m/s^2 (3 is default, might have to reduce this if the drone is destabilizing)

% quadrotor FRS and tracking error table
frs_filename = 'quadrotor_FRS_v7_a10_dt0.02.mat' ;
tbl_filename = 'quadrotor_tracking_error_table_dt0.02_vmax_5.25_zonotope.mat' ;
tracking_error_type = 'none' ; % 'table' or 'constant' or 'none' ;
tracking_error_constant_value = 0.1 ; %

% display
agent_camera_distance = 5 ; % default is 3
agent_camera_position = [-5;0;2] ; % default is [-3;0;1.5]
agent_camera_view = 'behind' ; % none, behind, above, or onboard
run_manual_mode_flag = false ;
verbose = 7 ;
plot_sim_flag = true ;

%% automated from here
% create world
W = zonotope_box_world_safe_rl('verbose',verbose,...
    'goal_type','box',...
    'use_wall_obstacles_flag',true) ;

% get agent, FRS, and tracking error table
A = quadrotor_agent('camera_view_style','behind','verbose',verbose,...
    'camera_view_style',agent_camera_view,...
    'camera_direction_offset',agent_camera_position,...
    'camera_follow_distance',agent_camera_distance) ;

disp('Loading quadrotor FRS')
FRS = load(frs_filename) ;

if strcmp(tracking_error_type,'table')
    if ~exist('tbl','var')
        disp('Loading quadrotor tracking error table')
        tbl = load(tbl_filename) ;
        tbl = tbl.tracking_error_table ;
    else
        disp('Table already loaded') ;
    end
else
    disp('Tracking error table does not need to be loaded') ;
    tbl = [] ;
end

% create agentHelper
AH = quadrotor_agentHelper(A,FRS,tbl,'v_max',v_max,'a_max',a_max,...
    'verbose',verbose,...
    'use_agent_for_initial_condition_flag',true,...
    'tracking_error_type',tracking_error_type,...
    'tracking_error_constant_value',tracking_error_constant_value) ;

% create rl simulator
S = rlsimulator(AH,W,'plot_sim_flag',plot_sim_flag,...
    'safety_layer','Z') ;

%%
S.reset() ;

% get the initial action
action_old = zeros(3,1) ;

% set the IsDone flag
IsDone = 0 ;

% run through keyboard presses
while (IsDone == 0) || (IsDone == 2)
    % convert delta_v into action
    action = 4*rand(3,1) - 2 ;
    action = bound_values(action,1) ;
    
    % step the rl simulator forward
    [~,~,IsDone,LoggedSignal,action_replaced,k_adjusted] = S.step(action) ;
    
    % check if the simulation is done
    if IsDone == 1 || IsDone == 3 || IsDone == 4 || IsDone == 5
        % 0. none of the below apply
        % 1. crash
        % 2. action replaced
        % 3. crash with safety layer on
        % 4. safely stopped but stuck
        % 5. reached goal!
        
        disp('-------------------')
        disp('SIMULATION COMPLETE')
        disp('-------------------')
        
        switch IsDone
            case 1
                disp('Crashed!')
            case 3
                disp('Crashed with safety layer on!')
            case 4
                disp('Safety stopped but stuck.')
            case 5
                disp('Reached goal!')
        end
        break
    end
end