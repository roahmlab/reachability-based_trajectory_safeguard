%% description
% This script tests the quadrotor_agentHelper.
%
% Note that you can use the keyboard to direct the quadrotor, by setting
% the run_manual_mode_flag to true. The controls are:
%   space --> no change in velocity
%   w     --> increase velocity forward
%   a     --> increase velocity left
%   s     --> increase velocity backward
%   d     --> increase velocity left
%   j     --> increase velocity up
%   k     --> increase velocity down
% In particular, each of these keys increases the quadrotor's commanded
% speed in the given direction by 1 m/s.
%wadw
% Author: Shreyas Kousik
% Created: 5 Aug 2020
% Updated: 7 Aug 2020
% clear;clc;close all;
plot_sim_flag = 1;
%% user parameters
% quadrotor FRS and tracking error table
frs_filename = 'quadrotor_FRS_v7_a10_dt0.02.mat' ;
tbl_filename = 'quadrotor_tracking_error_table_dt0.02_vmax_5.25_zonotope.mat' ;

% number of iterations to test
N_testing_iterations = 100 ;

% display
agent_camera_distance = 3%30 ; % default is 3
agent_camera_position =[-3;0;1.5]%[80;0;10] ; % default is [-3;0;1.5][-3;0;1.5]% 
agent_camera_view = 'behind' ; % none, behind, above, or onboard
run_manual_mode_flag = 0;
verbose = 0 ;

RTD_HLP_buffer = 0.1;
RTD_HLP_grow_tree_mode = 'keep';
%% automated from here
HLP = straight_line_HLP() ;
% create world
W = zonotope_box_world_safe_rl('verbose',verbose,...
    'goal_type','box',...
    'use_wall_obstacles_flag',true,...'N_tall',5,'N_wide',8,'N_long',0,'N_boxy',0) ;
'N_tall',20,'N_wide',20,'N_long',0,'N_boxy',0);
% get agent, FRS, and tracking error table
A = quadrotor_agent('verbose',verbose,...
    'camera_view_style',agent_camera_view,...
    'camera_direction_offset',agent_camera_position,...
    'camera_follow_distance',agent_camera_distance,'sensor_radius',10) ;
disp('Loading quadrotor FRS')
FRS = load(frs_filename) ;

if ~exist('tbl','var')
    disp('Loading quadrotor tracking error table')
    tbl = load(tbl_filename) ;
    tbl = tbl.tracking_error_table ;
else
    disp('table already loaded') ;
end

% create agentHelper
AH = quadrotor_agentHelper(A,FRS,tbl, HLP,'v_max',5,'a_max',3,'verbose',verbose) ;

%to fix the bug regarding vwelocity out of range
AH.use_agent_for_initial_condition_flag = true;
AH.tracking_error_type = 'constant';
A.move_method = 'integrator';

% create rl simulator
S = rlsimulator(AH,W,'plot_sim_flag',plot_sim_flag,...
    'safety_layer','Z') ;
AH.S = S;

if ~isempty(HLP)
AH.HLP.goal = S.W.goal;
AH.HLP.dimension = 3 ;
aginfo = AH.get_agent_info;
AH.HLP.setup(aginfo,S.W.get_world_info(aginfo))
end
%%
if run_manual_mode_flag
    % reset the simulator
    S.reset();
    
    % get the initial velocity
    v_cur = A.state(A.velocity_indices,end) ;
    
    % run through keyboard presses
    for idx = 1:N_testing_iterations
        
        % get keyboard press value
        keyboard_press_flag = waitforbuttonpress ;
        
        if keyboard_press_flag % 1 if keyboard was pressed
            value = double(get(gcf,'CurrentCharacter')) ;
            
            % create quadrotor velocity change based on keyboard press
            switch value
                case 32 % space --> no change
                    delta_v = [0;0;0] ;
                case 119 % w --> forward
                    delta_v = [1;0;0] ;
                case 97 % a --> left
                    delta_v = [0;1;0] ;
                case 115 % s --> backward
                    delta_v = [-1;0;0] ;
                case 100 % d --> right
                    delta_v = [0;-1;0] ;
                case 106 % j --> up
                    delta_v = [0;0;1] ;
                case 107 % k --> down
                    delta_v = [0;0;-1] ;
                otherwise
                    delta_v = [0;0;0] ;
            end
            
           
            
            % convert delta_v into action
            action = delta_v;
            
            % step the rl simulator forward
            [~,~,IsDone,LoggedSignal,action_replaced,k_adjusted] = S.step(action) ;
            
            % check if the simulation is done
            if IsDone == 1 || IsDone == 3 || IsDone == 4 || IsDone == 5
                % 1. crash       
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
            
            % reset v_cur to the adjusted parameter
            v_cur = k_adjusted ;
        end
    end
end