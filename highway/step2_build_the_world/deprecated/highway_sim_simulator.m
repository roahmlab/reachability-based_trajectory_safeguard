classdef highway_sim_simulator < handle
    % Class: simulator
    %
    % S = simulator(agents, worlds, planners, varargin)
    
    %% properties
    properties (Access = public)
        % basic properties
        agents = {agent()} ;
        worlds = {world()} ;
        planners = {planner()} ;
        
        N_agents = 1 ;
        N_worlds = 1 ;
        N_planners = 1 ;
        
        % user-friendly properties
        verbose = 1 ;
        
        % simulation
        max_sim_time = 100 ; % s per planner
        max_sim_iterations = 100 ; % per planner
        stop_count = 0 ;
        stop_threshold = 20 ;
        stop_sim_when_crashed = true ;
        allow_replan_errors = false ;
        save_planner_info = false ;
        manual_iteration = false ;
        collision_check_full_traj_after_sim_flag = false ;
        simulation_summary = [] ;
        
        
        % plotting
        figure_number = 1 ;
        figure_handle = [] ;
        plot_while_running = false ;
        debug_figure = 0;
        plotting_pause_time = 0.1 ; % s
        plot_order = 'WAP' ;
        save_gif = false ;
        save_gif_filename = 'simulator_gif_output.gif' ;
        save_gif_delay_time = 0.1 ;
        start_gif = true ;
        manually_resize_gif = true ;
        animation_time_discretization = 0.1 ; % s
        animation_linewidths = 1 ;
        clear_plot_before_animating_flag = false ;
        set_plot_linewidths_flag = false ;
        set_axes_while_animating_flag = false ;
        
        %RL needs stuff
        icur = 1;
        epscur = 0;
        W
        P
        A
        plot_in_loop_flag
        planning_time_vec
        stop_check_vec
        collision_check
        goal_check
        num_cars = 20
        small_num_cars_iteration = 1000000
        lane_width = 4 %m
        reset_counter = 0
        max_num_cars = 26
        min_num_cars = 10
        SIM_MAX_DISTANCE = 1000
        close_by = 0.01
        car_safe_dist = 26%40 % car 6 m long
        car_min_spd = 0
        car_max_spd = 5
        car_max_vis_dist = 120
        idm_dist_0 = 15
        idm_dist_1 = 30
        idm_dist_2 = 45
        idm_dist_3 = 60
        idm_dist_4 = 120
        
        no_coli_t = 3
        
        %action = 1 2 3 4 5  6 7 8 9  10 11 12;
        %         L L L L M  M M M R   R  R  R
        %        2d d 0 a 2d d m a 2d  d  m  a
        envCars % col 1 2  3 4    5    6
        %     x vx y vmax lane act
        
        A_MM = 0
        A_AM = 1
        A_BM = 2
        A_HM = 3
        A_ML = 4
        A_AL = 5% Not used at all
        A_BL = 6
        A_HL = 7 %not used often (<10)
        A_MR = 8
        A_AR = 9% not used at all
        A_BR = 10
        A_HR = 11 %not used often (<10)
        
        lev = '2';
        pref_act = -1;
        in_lane_thresh = 0.1;
        envCars_history = cell(200,1);
        inLane = true;
        U =[0;0];
        T= [0];
        T_brake = [];
        U_brake = [];
        Z_brake = [];

        D_mode = false;
    end
    
    %% methods
    methods
        %% constructor
        function S = highway_sim_simulator(agents, worlds, planners, varargin)
            % Constructor function: simulator
            %
            % Usage: S = simulator(agents, world, planners, varargin)
            %
            % This constructor takes in agents, which obey some physical
            % dynamics; world objects, which contains obstacles and goals;
            % and planners for the agent to perform receding-horizon traj.
            % planning in the worlds.
            
            % parse inputs
            S = parse_args(S,varargin{:}) ;
            
            % get number of agents, worlds, and planners
            S.N_agents = length(agents) ;
            S.N_worlds = length(worlds) ;
            S.N_planners = length(planners) ;
            
            % error if N_agents > 1 and N_agents ~= N_planners
            if S.N_agents > 1 && S.N_agents ~= S.N_planners
                error('Please provide either one agent, or one agent per planner')
            end
            
            % if the agent, world, or planner are alone, wrap them in cells
            if S.N_agents == 1 && ~iscell(agents)
                agents = {agents} ;
            end
            
            if S.N_worlds == 1 && ~iscell(worlds)
                worlds = {worlds} ;
            end
            
            if S.N_planners == 1 && ~iscell(planners)
                planners = {planners} ;
            end
            
            % wrap up construction
            S.agents = agents ;
            S.worlds = worlds ;
            S.planners = planners ;
        end
        function [newAct] = rplAct (S,eg_v, cf_d, cf_v, actAccl)
            if all(cf_d <= S.idm_dist_1) && all((eg_v-cf_v) >= 2)
                newAct = S.A_HM;
            elseif eg_v > S.car_min_spd && (all(cf_d <= S.idm_dist_0) || S.chkSafeDist(eg_v,cf_d,cf_v,actAccl)==false || (all(cf_d < S.idm_dist_2) && all((eg_v - cf_v) >= 2)))
                newAct =  S.A_BM;
            else
                newAct = S.A_MM;
            end
        end
        function [newAct] = rplAct_simple (S,eg_v, cf_d, cf_v, actAccl)
            if cf_d <= S.idm_dist_1 && (eg_v-cf_v) >= 2
                newAct = S.A_HM;
            elseif eg_v > S.car_min_spd && (cf_d <= S.idm_dist_0 || S.chkSafeDist_idxed(eg_v,cf_d,cf_v,actAccl,1)==false || (cf_d < S.idm_dist_2 && (eg_v - cf_v) >= 2))
                newAct =  S.A_BM;
            else
                newAct = S.A_MM;
            end
        end
        function [newAct] = rplAct_rear_simple (S,eg_v, cr_d, cr_v)
            if abs(cr_d)<= S.idm_dist_1 && (cr_v - eg_v) >= 2
                newAct = S.A_AM;
            else
                newAct = S.A_MM;
            end
            
        end
        function [newAct,reset_preference]= safeAct2(S,act,actAccl,inLane)
            eg_v = S.envCars(1,2);
            reset_preference = true;
            newAct = act;
            [forward_dist, forward_rel_spd] = S.getCar(4);%get ego car front info
            forward_spd = forward_rel_spd +eg_v;
            %             ego_lane = S.determine_ego_lane(2);
            %switch_lane_flg = 0;
            if inLane
                if any(act == [4 5 6 7])  && S.A.lane < 2 % change lane left
                    S.A.lane_des = S.A.lane + 1;
                elseif any(act ==[8 9 10 11]) && S.A.lane > 0% change lane right
                    S.A.lane_des = S.A.lane - 1;
                elseif act >= S.A_MM && act <= S.A_HM
                    S.A.lane_des = S.A.lane;
                else
                    newAct = S.rplAct(eg_v,forward_dist,forward_spd,actAccl);
                    return
                end
                
            end
            %                 ego_desired_lane = S.A.lane_des;
            
            % here this section helps with lane change, !Determine Ego lane should take into account of both
            if inLane
                if (~S.chkSafeDist_idxed(eg_v,forward_dist,forward_spd,actAccl,1))
                    newAct = S.rplAct(eg_v,forward_dist,forward_spd,actAccl);
                    return
                elseif forward_dist(1) >= S.car_max_vis_dist
                    if S.envCars(1,2)>31
                        newAct = S.A_MM;
                    else
                        newAct = S.A_AM;
                    end
                    return
                elseif S.A.lane == 0
                    if act >= S.A_MR
                        newAct = S.rplAct(eg_v,forward_dist,forward_spd,actAccl);
                        return
                    end
                elseif S.A.lane == 2
                    if act >= S.A_ML && act <= S.A_HL
                        newAct = S.rplAct(eg_v,forward_dist,forward_spd,actAccl);
                        return
                    end
                end
            end
            %%%%%% TODO: This following if can be simplified with 6-(S.A.lane*2+1)
            if S.A.lane_des == 0% change lane left
                %1   2
                %3 e 4
                %5   6
                %                 switch_lane_flg = S.A_ML;
                %                 switch_lane_rev = S.A_MR;
                [fdes_d,fdes_rel_v] = S.getCar_abs(6);
                fdes_v = fdes_rel_v + eg_v;
                [rdes_d,rdes_rel_v] = S.getCar_abs(5);
                rdes_v = rdes_rel_v + eg_v;
                %                 if inLane
            elseif  S.A.lane_des == 1%act >= S.A_MR && act <= S.A_HR% change lane right
                %1   2
                %3 e 4
                %5   6
                %                 switch_lane_flg = S.A_MR;
                %                 switch_lane_rev = S.A_ML;
                [fdes_d,fdes_rel_v] = S.getCar_abs(4);
                fdes_v = fdes_rel_v + eg_v;
                [rdes_d,rdes_rel_v] = S.getCar_abs(3);
                rdes_v = rdes_rel_v + eg_v;
            elseif S.A.lane_des == 2
                [fdes_d,fdes_rel_v] = S.getCar_abs(2);
                fdes_v = fdes_rel_v + eg_v;
                [rdes_d,rdes_rel_v] = S.getCar_abs(1);
                rdes_v = rdes_rel_v + eg_v;
            else
                error("lane des error")
            end
            
            if S.A.lane_des ~= S.A.lane
                switch_lane_flg = (S.A.lane_des - S.A.lane);
                if switch_lane_flg >= 1
                    switch_lane_flg = S.A_ML;
                    switch_lane_rev = S.A_MR;
                elseif switch_lane_flg <= 1 && switch_lane_flg >0
                    switch_lane_flg = S.A_MR;
                    switch_lane_rev = S.A_ML;
                    
                else
                    switch_lane_flg = 0;
                    switch_lane_rev = 0;
                    
                end
                
                [bk_d, bk_rel_v] = S.getCar_abs(6-(S.A.lane*2+1));%get ego car front info
                bk_v = bk_rel_v +eg_v;
                fcur_col=~S.chkSafeDist_idxed(eg_v, forward_dist(1), forward_spd(1),actAccl,1);
                fdes_col=~S.chkSafeDist_idxed(eg_v, fdes_d(1), fdes_v(1),actAccl,1);
                rcur_col =~S.chkRearSafe_idxed(eg_v,bk_d(1),bk_v(1),1);
                rdes_col=~S.chkRearSafe_idxed(eg_v,rdes_d(1),rdes_v(1),1);
                
                if S.plot_in_loop_flag
                    draw_offset = 8;
                    lane = S.determine_ego_lane(2);
                    %                      if fcur_col
                    plot(S.envCars(1,1) + draw_offset ,lane*4+2,'-s','MarkerSize',10,...
                        'MarkerEdgeColor','green',...
                        'MarkerFaceColor',[1 0 0.07]*fcur_col + [1  1 0]*~fcur_col)
                    %                      end
                    %                      if rcur_col
                    plot(S.envCars(1,1) - draw_offset ,lane*4+2,'-s','MarkerSize',10,...
                        'MarkerEdgeColor','green',...
                        'MarkerFaceColor',[1 0 0.07]*rcur_col + [1  1 0]*~rcur_col)
                    %                      end
                    %                      if switch_lane_flg == S.A_MR
                    %                         lane_add = -1;
                    %                      elseif switch_lane_flg == S.A_ML
                    %                          lane_add = 1;
                    %                      else
                    %                          lane_add = 0
                    %                      end
                    %
                    %                      if rdes_col
                    plot(S.envCars(1,1) - draw_offset ,(S.A.lane_des)*4+2,'-s','MarkerSize',10,...
                        'MarkerEdgeColor','red',...
                        'MarkerFaceColor',[1 0 0.07]*rdes_col + [1  1 0]*~rdes_col)
                    %                      end
                    %                       if fdes_col
                    plot(S.envCars(1,1) + draw_offset ,(S.A.lane_des)*4+2,'-s','MarkerSize',10,...
                        'MarkerEdgeColor','red',...
                        'MarkerFaceColor',[1 0 0.07]*fdes_col + [1  1 0]*~fdes_col)
                    %                      end
                    
                    
                end
                
                
                %                 fcur_v=forward_spd(1);
                %                 fdes_v=fdes_v(1);
                %                 fcur_d = forward_dist(1);
                %                 fdes_d = fdes_d(1);
                %                 rcur_d = bk_d(1);
                %                 rdes_d = rdes_d(1);
                %                 rcur_v = bk_v(1);
                %                 rdes_v = rdes_v(1);
                
                if inLane
                    if fdes_col || rdes_col || fcur_col
                        newAct = S.rplAct_simple(eg_v,forward_dist(1),forward_spd(1),actAccl);
                    else
                        newAct = act;
                        
                        %                         reference = false;
                    end
                else
                    
                    %                     if fdes_col && ~ rdes_col
                    %                         newAct = switch_lane_flg + S.rplAct_simple(eg_v,fdes_d,fdes_v,actAccl);
                    %                     elseif ~fdes_col && rdes_col
                    %                         newAct = switch_lane_flg + S.rplAct_rear_simple(eg_v,rdes_d,rdes_v);
                    %                     elseif ~fdes_col && ~ rdes_col && fcur_col
                    %                         newAct = switch_lane_flg + S.rplAct_simple(eg_v,fcur_d,fcur_v,actAccl);
                    %                     elseif fdes_col && rdes_col
                    %                         S.A.lane_des = S.A.lane;
                    %                         if fcur_col && ~ rcur_col
                    %                             newAct = switch_lane_rev + S.rplAct_simple(eg_v, fcur_d, fcur_v, actAccl);
                    %                         elseif ~ fcur_col && rcur_col
                    %                             newAct = switch_lane_rev + S.rplAct_rear_simple(eg_v, rcur_d, rcur_v);
                    %                         elseif ~ fcur_col && ~ rcur_col
                    %                             newAct = switch_lane_rev;
                    %                         else
                    %                             newAct = randi([4,11]);
                    %                         end
                    %                     else
                    %                         newAct = act;
                    %                         reset_preference = false;
                    %                     end
                    if fdes_col || rdes_col
                        S.A.lane_des = S.A.lane;
                        newAct = switch_lane_rev;
                    else
                        newAct = act;
                        reset_preference = false;
                    end
                    
                end
                
            end
            
            %             if switch_lane_flg&&(~S.chkSafeDist(eg_v, fdes_d, fdes_v,actAccl) || ~S.chkRearSafe(eg_v,rdes_d,rdes_v))
            %                     newAct = S.rplAct(eg_v,forward_dist,forward_spd,actAccl);
            %             end
        end
        
        function [Safe] = chkSafeDist(S,ego_v,f_d_a, f_v_a, actAccl)
            safe_a = zeros(2,1);
            for i = 1: 2
                f_d = f_d_a(i);
                f_v = f_v_a(i);
                max_ego_spd = min([S.car_max_spd,ego_v+actAccl*2]);%bool*a
                min_front_spd = max([S.car_min_spd,f_v-4]);
                
                if max_ego_spd > min_front_spd
                    Safe = (f_d - 15) - 4*(max_ego_spd - min_front_spd) > 0; %15m of dist at least, 4 seconds of time to collision
                    %                     if ~Safe
                    %                         warning("Not enough dist in front");
                    %                     end
                else
                    Safe = true;
                end
                safe_a(i) = Safe;
            end
            Safe = all(safe_a);
        end
        function [Safe] = chkSafeDist_idxed(S,ego_v,f_d_a, f_v_a, actAccl,idx)
            i =  idx;
            f_d = f_d_a(i);
            f_v = f_v_a(i);
            max_ego_spd = min([S.car_max_spd,ego_v+actAccl*2]);%bool*a
            min_front_spd = max([S.car_min_spd,f_v-4]);
            
            if max_ego_spd > min_front_spd - 1
                Safe =( (f_d - 15) - 3*(max_ego_spd - min_front_spd) )> 0; %15m of dist at least, 3 seconds of time to collision
            else
                Safe = true;
            end
        end
        function [Safe] = chkRearSafe(S,ego_v,r_d_a, r_v_a)
            safe_a = zeros(2,1);
            for i = 1: 2
                r_d = r_d_a(i);
                r_v = r_v_a(i);
                min_ego_spd = max([S.car_min_spd,ego_v-4]);
                
                if min_ego_spd < r_v
                    Safe = ((r_d+15) + (r_v - min_ego_spd) * 2 ) < 0;
                elseif min_ego_spd == S.car_min_spd && r_v == S.car_min_spd
                    Safe = (r_d+10 ) < 0;
                else
                    Safe = r_d + 15 <0;
                end
                safe_a(i) = Safe;
            end
            Safe = all(safe_a);
        end
        function [Safe] = chkRearSafe_idxed(S,ego_v,r_d_a, r_v_a,idx)
            i =idx;
            r_d = r_d_a(i);
            r_v = r_v_a(i);
            min_ego_spd = max([S.car_min_spd,ego_v-4]);
            
            if min_ego_spd < r_v
                Safe = ((r_d+15) + (r_v - min_ego_spd) * 2 ) < 0;
            elseif min_ego_spd == S.car_min_spd && r_v == S.car_min_spd
                Safe = (r_d+10 ) < 0;
            else
                Safe = r_d + 15 <0;
            end
            
        end
        %% run simulation
%         function examine_frs(S)
            %reset enviroment with left/right extreme initial conditions
            %plot FRS 0-0.5 in red, 0.5-3 in green
            %move agents
            
            %plot agent
            
            
%             for h_idx = 1: hr
%                 for w_idx = 4:6
%                     for v_idx = vr
%                         for vd_idx = 1:num_v_des
%                             S.icur=1;
%                             S.envCars_history = cell(200,1);
%                             iniOb = zeros(21,1);
%                             LoggedSignals= struct;
%                             %         t_max = S.max_sim_time ;
%                             iter_max = S.max_sim_iterations ;
%                             S.plot_in_loop_flag = S.plot_while_running;
%                             S.W = S.get_world(widx) ;
%                             % get agent and planner
%                             S.A = S.get_agent(pidx) ;
%                             S.P = S.get_planner(pidx) ;
%                             
%                             % get agent and world ready
%                             S.W.reset() ;
%                             
%                             %randomize start location for ego
%                             start_lane = 1;
%                             y_ini = 4*start_lane + 2;
%                             S.A.lane = start_lane;
%                             S.A.lane_des = start_lane;
%                             start_v = 22;
%                             if isprop(S.A,'desired_initial_condition')
%                                 S.A.desired_initial_condition(2) = y_ini;
%                                 S.A.desired_initial_condition(4) = start_v;
%                                 S.A.reset([S.A.desired_initial_condition]) ;
%                             else
%                                 S.A.reset(S.W.start) ;
%                             end
%                             
%                             
%                             
%                             % reset other cars on the road
%                             S.reset_counter = S.reset_counter + 1;
%                             if S.reset_counter > S.small_num_cars_iteration
%                                 S.num_cars = S.min_num_cars + randi(S.max_num_cars-S.min_num_cars);
%                             end
%                             
%                             S.placeCars()
%                             
%                             S.W.setup(S.envCars)
%                             
%                             % get planner ready
%                             agent_info = S.A.get_agent_info() ;
%                             world_info = S.W.get_world_info(agent_info);
%                             agent_info.lane = start_lane;
%                             agent_info.spd = start_v;
%                             
%                             S.P.setup(agent_info,world_info) ;
%                             
%                             % check to make sure gif start is ready
%                             if S.save_gif
%                                 S.start_gif = true ;
%                             end
%                             
%                             % initialize plot
%                             %             if  S.plot_while_running
%                             %                 S.plot(widx,pidx)
%                             %             end
%                             
%                             % preallocate for storing planning time spent
%                             S.planning_time_vec = nan(1,iter_max) ;
%                             
%                             % reset the stop counter
%                             S.stop_count = 0 ;
%                             S.stop_check_vec = false(1,iter_max) ;
%                             
%                             % reset the crash and goal checks just in case
%                             S.collision_check = false ;
%                             S.goal_check = false ;
%                             
%                             %             iniOb = agent_info.state(:,1);
%                             LoggedSignals.episode_result = 'initial';             
%                             S.step(4,1,1);
%                         end
%                     end
%                 end
%             end  
%         end
        
        function [Observation,Reward,IsDone,LoggedSignals]= step(S,action,widx,pidx)
            % observation: distance to goal, rel_position(x, y),
            % position (x, y)_ego, current_velocity, average_velocity
            IsDone = 0;
            Observation = zeros(5,1);   %%%%%change
            
            LoggedSignals = struct;
            
            %% get agent info
            agent_info = S.A.get_agent_info() ;
            
            %% get world info
            % given the current state of the agent, query the world
            % to get the surrounding obstacles
            world_info = S.W.get_world_info(agent_info,S.P);
            %% replan
            % given the current state and obstacles, query the
            % current planner to get a control input
            t_plan_spent = tic ;
            newAct = action;       %% new action  2 -4
            
            if S.D_mode == true
                if any(newAct == [0 4 8])%M
                    agent_info.vx_des = 0;
                elseif any(newAct == [1])%A
                    agent_info.vx_des = 2;
                elseif any(newAct == [2])%B
                    agent_info.vx_des = -2;
                elseif any(newAct == [3])%H
                    agent_info.vx_des = -4;
                else
                end
                
                if any(newAct == [4])
                    agent_info.y_des = 1;
                elseif any(newAct == [8])
                    agent_info.y_des = -1;
                else
                    agent_info.y_des = 0;
                end
            else
                newAct(1) = (action(1)+1)/2*6-4;   %%%%%%  vx
                newAct(2) = action(2)*1;%+1)/2*8+2;    %%%%%%  y
                agent_info.vx_des = newAct(1);          % -4, 2
                agent_info.y_des = newAct(2);      %-4,4
            end
            %             LoggedSignals.agent_info =agent_info;
            %             LoggedSignals.world_info =world_info;
            %             LoggedSignals.envCars = S.envCars;
            
            %             end
            actAccl = 0;
            %             if newAct == S.A_AM || newAct == S.A_AL || newAct == S.A_AR
            %                 actAccl = 1;
            %             else
            %                 actAccl = 0;
            %             end
            
            agent_info.lane = S.A.lane;
    
            
       
            if  S.P.mode == 'N'||S.P.mode == 'Z'      %change   N no safety
                reset_prefer = true;
                if S.D_mode == false
                    agent_state = agent_info.state(:,end);
                end

                if agent_info.vx_des > 2
                    agent_info.vx_des = 2;
                elseif agent_info.vx_des < -4
                    agent_info.vx_des = -4;
                end

                if agent_info.y_des > 1
                    agent_info.y_des = 1;
                elseif agent_info.y_des < -1
                    agent_info.y_des = -1;
                end
            end
            
            if S.P.mode == 'R'
                if newAct ~= action % here it may be ok to use safeact and still have some reward, look at paper
                    Reward = -1;
                    IsDone = 2;
                    if S.plot_in_loop_flag
                        %figure(1)
                        %rectangle('Position', [S.envCars(1,1)-100, 0, 300, 13], 'EdgeColor', [1, 0, 0, 0], 'FaceColor', [1, 0, 0, 0.1], 'LineWidth', 1);
                        %rectangle('Position', [S.envCars(1,1)-100, 0, 300, 13], 'EdgeColor', [1, 0, 0, 0], 'FaceColor', [1, 0, 0, 0.1], 'LineWidth', 1);
                        
                    end
                    
                    LoggedSignals.replaced_action = 'true';
                    LoggedSignals.secondexp = newAct;
                    
                else
                    LoggedSignals.secondexp = NaN;
                end
            end
           
            
            if S.allow_replan_errors         % safety layer %%%%%%%%%%%
                [T_nom,U_nom,Z_nom,replaced,K] = S.P.replan(agent_info,world_info) ; % estimate next movement T position,  U norm
                %%%%%%%%%%
                if ~isempty(T_nom)
                    S.T_brake = T_nom(101:end)-T_nom(101);
                    S.U_brake = U_nom(:,101:end);
                    S.Z_brake = Z_nom(:,101:end);
                end
                if replaced
                    IsDone = 2;
                    if isempty(K)
                            LoggedSignals.replaced_action = 'not_found';
                            LoggedSignals.secondexp = single([0; 0]);
                        else
                            LoggedSignals.replaced_action = 'true';
                            LoggedSignals.secondexp = single(reshape(K,[2,1]));
                        end
                    else
                        LoggedSignals.replaced_action = 'false';
                        LoggedSignals.secondexp = single([0; 0]);
                end
            else
                try
                    [T_nom,U_nom,Z_nom,~,~] = S.P.replan(agent_info,world_info) ;
                catch
                    S.vdisp(['Planner ',num2str(pidx),' errored while ',...
                        'replanning!'])
                    T_nom = [] ; U_nom = [] ; Z_nom = [] ;
                end
            end
            
            t_plan_spent = toc(t_plan_spent) ;
            %         planning_time_vec(icur) = t_plan_spent ;
            S.vdisp(['Planning time: ',num2str(t_plan_spent),' s'],4)
            
            if S.P.mode == 'Z'
                reset_prefer =S.P.reset_prefer;
                %                 if reset_prefer
                %                     Reward = -1;
                %                     IsDone = 2;
                %                     LoggedSignals.replaced_action = 'true';
                %                     LoggedSignals.secondexp = S.P.class_num;
                %                     newAct = S.P.class_num;
                %                 else
                %                     LoggedSignals.secondexp = NaN;
                %                 end
            end
            
            
            %plot after getting plans since world info has info about
            %env car predictions.
            %             if S.plot_in_loop_flag
            %                 %S.plot(widx,pidx)
            %
            %                 %                 scatter(S.envCars(:,1),S.envCars(:,3))
            %                 if S.save_gif
            %                     error('Shreyas this is unfinished!')
            %                 else
            %                     pause(S.plotting_pause_time) ;
            %                 end
            %             end
            %%
            %             lane_ego=S.determine_ego_lane(1);
            %             desired_lane_ego = lane_ego+(S.P.y_des- agent_info.state(2,end) )%(S.P.y_des - 2)/4;
            
            %% move agent
            % update the agent using the current control input, so
            % either stop if no control was returned, or move the
            % agent if a valid input and time vector were returned
            if size(T_nom,2) < 2 || size(U_nom,2) < 2 || T_nom(end) == 0
                t_move = 4;
                S.vdisp('cannot find replacement action, stopping!',2)
                %S.setAction(); %%%%%%%%%%%%%%%%%%%%%  move enviroment car
                %S.applyAction();
                S.A.move(t_move,S.T_brake,S.U_brake,S.Z_brake) ; %%%%%%%% move ego car
                S.envCars(1,:)= [S.A.state(1,end), S.A.state(4,end), S.A.state(2,end), 32,S.A.lane_des ,S.A_MM]; %%%%%%%%%%%% x , vx , y , 32 ,  lane , action       ego car parameter
                S.W.update_envCars(S.envCars) ;%this is for updating states of other obstacles for plotting and collision check
                S.envCars_history{S.icur} = S.envCars;
                S.icur = S.icur +1;
                S.T_brake= [0 10];
                S.U_brake = [0 0; S.U_brake(2,end) S.U_brake(2,end)];
                S.stop_count = S.stop_count + 1;
                if S.stop_count > 3
                    IsDone = 4;
                end
            else
                t_move = 2;
                %S.setAction(); %%%%%%%%%%%%%%%%%%%%%  move enviroment car
                %S.applyAction();
                S.A.move(t_move,T_nom,U_nom,Z_nom) ; %%%%%%%% move ego car
                S.envCars(1,:)= [S.A.state(1,end), S.A.state(4,end), S.A.state(2,end), 32,S.A.lane_des ,S.A_MM]; %%%%%%%%%%%% x , vx , y , 32 ,  lane , action       ego car parameter
                S.W.update_envCars(S.envCars) ;%this is for updating states of other obstacles for plotting and collision check
                S.envCars_history{S.icur} = S.envCars;
                S.icur = S.icur +1;
            end 
                min_index = 2;
                min_dis = inf;
                min_y = 0;
                
                for idx_neighbor = [2,4,6] % observation here can do better, use a hybrid of desired pose and actual pose;
                    %currently using desired pose since it will eventually converge to desired pose and current pose doesnt really matter
                    [rel_dist,rel_v,rel_y, ego_vy]=S.W.getCar_abs(idx_neighbor,agent_info);       %%%%%%%%%%%%%%%%%%%%%%%%% get car get_abs
                    if rel_dist < min_dis
                        min_index = idx_neighbor;
                        min_dis = rel_dist;
                        min_y = rel_y;
                    end
                    

%                     Observation(3*(idx_neighbor-1)+1:3*idx_neighbor) = [rel_dist;rel_v;rel_y];
                end
                Observation(1) = min_dis;
%                 if min_y >= 0
%                     Observation(2) = 1/(min_y+1);
%                 else 
%                     Observation(2) = 1/(min_y-1);
%                 end
                Observation(2) = min_y;
                %Observation(19:21) = [S.envCars(1,2);S.envCars(1,3);ego_vy];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                Observation(3) = S.envCars(1,1); %ego_x
                Observation(4) = S.envCars(1,3); %ego_y
                Observation(5) = S.envCars(1,2); %ego_vx
                %                 scaled_Ob = S.get_scaled_indicator(Observation);
                
                [thresh, lane] = min(abs([2, 6, 10]- S.P.y_des));
                lane = lane-1;
                if S.inLane
                    if  thresh>0.3% in lane
                        S.inLane = false;
                        
                    end
                else
                    if thresh< 0.2
                        S.inLane = true;
                        S.A.lane = lane;
                    end
                end
                
                if S.plot_in_loop_flag
                    if S.debug_figure
                        figure(5);clf; hold on;
                        
                        for idx_neighbor = 1:6 % observation here can do better, use a hybrid of desired pose and actual pose;
                            %currently using desired pose since it will eventually converge to desired pose and current pose doesnt really matter
                            scatter(Observation(3*idx_neighbor-2),Observation(3*idx_neighbor));
                            text(Observation(3*idx_neighbor-2)+1,Observation(3*idx_neighbor),num2str(Observation(3*idx_neighbor-1),3))
                            xline(S.car_max_vis_dist,'--')
                            xline(-S.car_max_vis_dist,'--')
                            xlim([-150 150])
                            
                            %                             Observation(3*(idx_neighbor-1)+1:3*idx_neighbor) = [rel_dist;rel_v;rel_y];
                        end
                        scatter(0,0)
                    end
                    figure(1)
                    S.plot(widx,pidx)
                    rectangle('Position', [S.envCars(1,1)-20,-5, 4, 4], 'EdgeColor', 'b', 'FaceColor', [1, 0, 0, 0], 'LineWidth', 1);
                    plot(S.envCars(1,1)-18+(S.P.vx_des-S.envCars(1,2)), -3+S.P.y_des-S.envCars(1,3),'-s','MarkerSize',10,...
                        'MarkerEdgeColor','blue',...
                        'MarkerFaceColor',[200/255*reset_prefer 0 1])
                    for i = 2:size(S.envCars,1)
                        text (S.envCars(i,1)-2,S.envCars(i,3), num2str(S.envCars(i,2),3));
                    end
                    %text (S.envCars(1,1)+3, S.P.y_des, num2str(S.P.vx_des));
                    xline(S.envCars(1,1)+S.car_max_vis_dist,'--')
                    xline(S.envCars(1,1)-S.car_max_vis_dist,'--')
                    text (S.envCars(1,1)+100,6, "N="+num2str(newAct,2));
                    text (S.envCars(1,1)+100,2, "A="+num2str(action,2))
                    text (S.envCars(1,1)+90,2, "D="+num2str(IsDone));
                    text (S.envCars(1,1)-100,6, "V="+num2str(agent_info.state(4,end)));
                    text (S.envCars(1,1)-100,2, "Y="+num2str(agent_info.state(2,end)));
                    text (S.envCars(1,1)-100,2-4, "Ydes="+num2str(S.P.y_des));
                    text (S.envCars(1,1)-100,6+4, "Vdes="+num2str(S.P.vx_des));
                    
                    %if S.debug_figure
%                     S.U = [S.U;U_nom(1,1)];
%                         figure(6) ;plot(diff(S.A.state(4,:))./diff(S.A.time));title('Acceleration');figure(1);
%                         figure(7) ;plot(S.U);title('Ref Velocity');
%                         figure(1);
                    %end
                    %                     text (S.envCars(1,1)+3, S.P.y_des, num2str(S.P.vx_des));
                    %                     text (S.envCars(1,1)+100,6, "N="+num2str(newAct));
                    %                     text (S.envCars(1,1)+100,2, "A="+num2str(action))
                    %                     text (S.envCars(1,1)+90,2, "D="+num2str(IsDone));
                    %                     text (S.envCars(1,1)-100,6, "V="+num2str(agent_info.state(4,end)));
                    %                     text (S.envCars(1,1)-100,2, "Y="+num2str(agent_info.state(2,end)));
                    % text (S.envCars(1,1)+90,2, "D="+num2str(IsDone));
                end
                
            
            if IsDone == 0
                Reward = S.getRew(Observation);% 1-sum((world_info.goal-agent_info.state(1:2,end)).^2)/sum((world_info.goal-world_info.start(1:2)).^2);%%%%%%%%%%%%
                
            elseif IsDone == 2        %%%%%%%%%%%%%%%%%%% get reward only
                Reward = S.getRew(Observation);  %%%%%LoggedSignals.newreward need to use
            else
                Reward = 0;
            end
            if S.plot_in_loop_flag
                text (S.envCars(1,1)+100,10, "R="+num2str(Reward, 2));
            end
            %cleanup
            %                 lane_num = lane_num -1 ;
            if reset_prefer || S.inLane
                S.pref_act = -1;
            else
                if newAct == 4 || newAct == 5 || newAct == 6 || newAct == 7
                    S.pref_act = 4;
                elseif  newAct == 8 || newAct == 9 || newAct == 10 || newAct == 11
                    S.pref_act = 8;
                end
            end
            
            
            
            
            %% Note (22 July 2019)
            % Dynamic obstacles are treated as follows:
            %   1) W.get_world_info should return a prediction
            %   2) the agent is moved according to the prediction
            %   3) the world moves the obstacles (according to the
            %      agent's movement data if needed) and then checks
            %      for collisions in W.collision_check
            
            %% crash and goal check
            % check if the agent is near the desired goal or if it
            % crashed
            S.vdisp('Checking if agent reached goal or crashed...',3)
            agent_info = S.A.get_agent_info() ;
            S.goal_check = S.W.goal_check(agent_info) ;
            S.collision_check = S.W.collision_check(agent_info,false) ;
            if isa(S.A,'multi_link_agent')
                S.vdisp('Checking for self-intersection.',2)
                S.collision_check = S.collision_check || S.A.self_intersection_flag ;
            end
            
            if S.goal_check
                Reward = Reward + 30 * (norm(S.W.goal-S.W.start(1:2)))/S.A.time(end);
                IsDone = 5;
            end
            % goal has no meaning here, just want cruising
            %             if S.goal_check
            %                 S.vdisp('Reached goal!',2) ;
            %                 IsDone = 1;
            %                 Reward = 1;
            %                 LoggedSignals.episode_result = 'goal reached';
            %             end
            
            if S.collision_check && S.stop_sim_when_crashed
                if S.plot_in_loop_flag
                    S.plot(widx,pidx)
                end
                S.vdisp('Crashed!',2) ;
                IsDone = 1;
                Reward = Reward;  %%%%%%%%%%%%%%%%
                LoggedSignals.episode_result = 'collision';
                if  (S.P.mode == 'Z' && reset_prefer)
                    IsDone = 3;
                    Reward = Reward;% discard this exp since new_act caused crash %%%%%%%%%%%%%%%%
                elseif S.P.mode == 'R' && newAct ~= action
                    IsDone = 3;
                    Reward = Reward;% discard this exp since new_act caused crash %%%%%%%%%%%%%%%%
                end
            else
                LoggedSignals.episode_result = 'safe';
            end
            
            
                
            
            if  IsDone == 1 || IsDone == 3 ||IsDone == 4     %%%%%%%%%%%%%%%% replay
                agent_name_cell{pidx} = S.A.name ;
                planner_name_cell{pidx} = S.P.name ;
                %                     trajectory_cell{pidx} = Z ;
                total_simulated_time_cell{pidx} = T_nom ;
                control_input_cell{pidx} = U_nom ;
                %                     control_input_time_cell{pidx} = TU ;
                %                     total_real_time_cell{pidx} = runtime ;
                total_iterations_cell{pidx} = S.icur-1;
                %                     planning_times_cell{pidx} = planning_time_vec ;
                collision_check_cell{pidx} = S.collision_check ;
                %                     goal_check_cell{pidx} = goal_check ;
                %                     stop_check_cell{pidx} = stop_check_vec ;
                t_plan_cell{pidx} = S.P.t_plan ;
                t_move_cell{pidx} = S.P.t_move ;
                agent_info_cell{pidx} = agent_info ;
                %                     obstacles_cell{pidx} = S.W.obstacles;
                
                
                %                                  'trajectory',trajectory_cell,...
                %                                  'total_real_time',total_real_time_cell,...
                %                                  'planning_time',planning_times_cell,...
                %                                  'goal_check',goal_check_cell,...
                %                                  'stop_check',stop_check_cell,...
                %                                  'control_input_time',control_input_time_cell,...
                %                                  'planner_info',planner_info_cell,...
                %                                  'obstacles',obstacles_cell,...
                %                                  'planner_indices',planner_indices,...
                %                                  'N_obstacles',S.W.N_obstacles,...
                %                                  't_max',t_max,...
                %                                  'iter_max',iter_max,...
                summary= struct('agent_name',agent_name_cell,...
                    'planner_name',planner_name_cell,...
                    'total_iterations',total_iterations_cell,...
                    'collision_check',collision_check_cell,...
                    'total_simulated_time',total_simulated_time_cell,...
                    'agent_info',agent_info_cell,...
                    't_plan',t_plan_cell,...
                    'control_input',control_input_cell,...
                    't_move',t_move_cell,...
                    'start',S.W.start,...
                    'goal',S.W.goal,...
                    'bounds',S.W.bounds,...
                    'notes','') ;
                LoggedSignals.summary = summary;
                LoggedSignals.envCars_hist = S.envCars_history;
                
                Filename = sprintf('sim_summary_%c_%s_%s-%s.mat', [S.P.mode IsDone S.epscur strcat(num2str(S.icur),"_",datestr(now,'HH-MM-SS.FFF'))]);
                fid = fopen(Filename,"wt")
                fclose(fid);
                %save(Filename,"LoggedSignals")
            end
            
           %% convert back to isdone = 0 or 1
%             if IsDone == 0
%                 LoggedSignals.episode_result = 'Executed Network Action';
%                 IsDone = 0;
%             elseif IsDone == 1
%                 LoggedSignals.episode_result = 'Crashed';
%                 IsDone = 1;
%             elseif IsDone == 2
%                 LoggedSignals.episode_result = 'Action Edited by Safety Layer';
%                 IsDone = 0;
%             elseif IsDone == 3
%                 LoggedSignals.episode_result = 'Action Edited and Crashed';
%                 IsDone = 1;
%             elseif IsDone == 4
%                 LoggedSignals.episode_result = 'Cannot find Replacement Action';
%                 IsDone = 1;
%             elseif IsDone == 5
%                 LoggedSignals.episode_result = 'Reached Goal';
%                 IsDone = 1;
%             end
            
        end
        function applyAction(S)
            t = S.P.t_move;
            S.envCars(:,1) = S.envCars(:,1) + S.envCars(:,2)*t;
            out_of_bounds_idx = abs(S.envCars(:,1) - S.A.state(1,end)) > S.SIM_MAX_DISTANCE;% x pos of ego
            S.envCars(out_of_bounds_idx,1) = S.envCars(out_of_bounds_idx,1) - sign(S.envCars(out_of_bounds_idx,1) - S.A.state(1,end))*2*S.SIM_MAX_DISTANCE;
            S.envCars(:,2) =S.envCars(:,2) + t*((S.envCars(:,6) == S.A_HM) * (-4) + (S.envCars(:,6) == S.A_BM) * (-2) + (S.envCars(:,6) == S.A_AM) * 2 );
            S.envCars(S.envCars(:,2) < 22,2) = 22;
            S.envCars(S.envCars(:,2) > S.envCars(:,4),2) = S.envCars(S.envCars(:,2) > S.envCars(:,4),4);
        end
        function [lane_ego] = determine_ego_lane(S,i)
            if i == 1
                ego_y = S.A.state(2,end);
            else
                lane_ego = S.A.lane;
                return
            end
            [~,lane_ego]=min([abs(ego_y-2),abs(ego_y-6),abs(ego_y-10)]);% 1,2,3
            lane_ego = lane_ego-1;%0,1,2
        end
        function placeCars(S)
            S.envCars = zeros(S.num_cars, 6);
            lane_ego = S.determine_ego_lane(2);
            spd_limit = floor(S.A.state(4,end)/2) *2;
            S.envCars(1,:) = [S.A.state(1,end), S.A.state(4,end), S.A.state(2,end), spd_limit, lane_ego ,S.A_MM];
            %envCars % col 1 2  3 4    5    6
            %     x vx y vmax lane act
            MAX_DIST =  S.SIM_MAX_DISTANCE;
            for i = 2:S.num_cars
                laneOverlap = true;
                while laneOverlap
                    xPos = MAX_DIST * (rand - 0.5) + S.envCars(1,1)+50;
                    laneNum = randi([0,2]);
                    laneOverlap = false;
                    for j = 1:i
                        %if overlap, y dir close and x dir close
                        % TODO: check with eco car as well, ego not here
                        if  abs(xPos - S.envCars(j,1)) < S.car_safe_dist %not useful any more, static obstacles%&&laneNum ==S.envCars(j,5)
                            laneOverlap = true;
                            break;
                        end
                    end
                    if ~laneOverlap % assign car to new slot
                        S.envCars(i,1) = xPos;
                        S.envCars(i,2) = 0;%S.car_min_spd+rand*(S.car_max_spd-S.car_min_spd);
                        S.envCars(i,3) = 4*laneNum + 2;
                        S.envCars(i,4) = 0;%S.car_max_spd-rand*(S.car_max_spd-S.car_min_spd)/3;
                        S.envCars(i,5) = laneNum;
                        S.envCars(i,6) = S.A_MM;%act #7 is maintain, for env cars only some actions are available
                        %action = 1 2 3 4 5  6 7 8 9  10 11 12;
                        %         L L L L M  M M M R   R  R  R
                        %        2d d 0 a 2d d m a 2d  d  m  a
                    end
                end
            end
        end
        function placeCars_blocked(S)
            S.envCars = zeros(S.num_cars, 6);
            lane_ego = S.determine_ego_lane(2);
            spd_limit = floor(S.A.state(4,end)/2) *2;
            S.envCars(1,:) = [S.A.state(1,end), S.A.state(4,end), S.A.state(2,end), spd_limit, lane_ego ,S.A_MM];
            %envCars % col 1 2  3 4    5    6
            %     x vx y vmax lane act
            MAX_DIST =  410;
            for i = 2:S.num_cars
                
                    xPos = -MAX_DIST*rand;
                    laneNum = randi([0,2]);
                    laneOverlap = false;
                    
                    if ~laneOverlap % assign car to new slot
                        S.envCars(i,1) = xPos;
                        S.envCars(i,2) = 0;%S.car_min_spd+rand*(S.car_max_spd-S.car_min_spd);
                        S.envCars(i,3) = 4*laneNum + 2;
                        S.envCars(i,4) = 0;%S.car_max_spd-rand*(S.car_max_spd-S.car_min_spd)/3;
                        S.envCars(i,5) = laneNum;
                        S.envCars(i,6) = S.A_MM;%act #7 is maintain, for env cars only some actions are available
                        %action = 1 2 3 4 5  6 7 8 9  10 11 12;
                        %         L L L L M  M M M R   R  R  R
                        %        2d d 0 a 2d d m a 2d  d  m  a
                    end
                
            end
        end
        function [Observation]=get_scaled_indicator(S,Observation)
            max_vy_ego = 1.5;
            for i = 1:6
                idx = (i-1)*4+1;
                Observation(idx)  = Observation(idx)/ S.car_max_vis_dist;
                idx = (i-1)*4+2;
                Observation(idx)  = Observation(idx)/ (S.car_max_spd-S.car_min_spd);
                idx = (i-1)*4+3;
                Observation(idx)  = Observation(idx)/ 8; % full road width
                idx = (i-1)*4+4;
                Observation(idx)  = Observation(idx)/ max_vy_ego; % assume y speed < 3m/s
            end
            Observation(25) = Observation(25) / (S.car_max_spd-S.car_min_spd); %lane 3 center
            % [S.envCars(1,2);S.envCars(1,3); -rel_vy];
            Observation(26) = Observation(26) / 12;
            Observation(27) =  Observation(27) /max_vy_ego;
        end
        %
        
        function [total_rew]= getRew(S,observation)
                %observation(2)
                y_near_obst = observation(2)+ observation(4);
                
                if (observation(4) <= 10) && (observation(4) >= 2)
                    Reward_Y = 5 * exp(- 1/(observation(4).^2+1));
                elseif (observation(4) < 2)
                    Reward_Y = 5 * exp(- 1/((y_near_obst-2)^2+1)) - 2 * abs(abs(observation(4))-2) ;
                elseif (observation(4) > 10)
                    Reward_Y = 5 * exp(- 1/((y_near_obst-10)^2+1)) - 2 * abs(abs(observation(4))-10) ;
                end
                %Reward_Y = 5 * exp(- 1/(observation(2)^2+1));
                %* exp(- 1/((observation(1)+1)^2 / 8))
                %Reward_Y
                %Reward_punish = 30 / ((norm(S.W.goal - [observation(4);observation(5)]) / 1000)+1);

%                 if (observation(5)==0)
%                         Reward_speed = -2;
%                 elseif (observation(5)<1.5)
%                         Reward_speed = -0.2 *observation(5);
%                 elseif (observation(5)>1.5) && (observation(5)<3)
%                         Reward_speed = 0.1*(observation(5));
%                 else 
%                         Reward_speed = 0.2 * observation(5); 
%                 end 
                Reward_speed = 4 * exp(- 1./(observation(5)^2+1)) - 2;
          
                total_rew = Reward_Y  + Reward_speed;  % Reward_punish
                
            
%             total_rew = S.envCars(1,2)/5+(1000-norm(S.W.goal - [S.envCars(1,1);S.envCars(1,3)]))/1000;
        
%             [dist, rel_v] = S.getCar(4); % fron car speed -ego
%             lane = S.determine_ego_lane(1);
%             dist= dist(1); rel_v = -rel_v(1);
%             safe_dist = 15 + rel_v * S.no_coli_t;
%             midDist =  safe_dist;
%             rewX = 1;
%             if dist < S.car_max_vis_dist - 1 && rel_v > 0
%                 if S.plot_in_loop_flag
%                     figure(1)
%                     plot_y = lane *4 ;
%                     plot_x = S.envCars(1,1) + dist - midDist;
%                     rectangle('Position', [plot_x, plot_y, midDist, 4], 'EdgeColor',[1, 0, 0, 0] , 'FaceColor', [1, 0, 0, 0.1], 'LineWidth', 1);
%                     
%                 end
%                 if  dist < midDist
%                     rewX = exp(-(dist-midDist)^2/(10*midDist));
%                     if S.plot_in_loop_flag
%                         x_ary = 0:0.5:midDist;
%                         y_ary = exp(-(x_ary-midDist).^2./(10*midDist));
%                         if S.debug_figure
%                             figure(4); clf;hold on;set(gca, 'XDir','reverse')
%                             
%                             plot(x_ary,y_ary);
%                             plot(dist, rewX,'r*');
%                             title('Xrew')
%                         end
%                     end
%                     
%                 end
%             else
%                 if S.plot_in_loop_flag
%                     if S.debug_figure
%                         figure(4); clf;
%                     end
%                     
%                 end
%             end
%             %%
%             distY  = lane*4+2; % default  y
%             if dist > S.car_max_vis_dist - 1
%                 maxV = 32;
%             else
%                 maxV = S.car_min_spd;
%                 %choose the best lane and speed of the car in front
%                 for i = 1:3
%                     [distf, rel_vf] = S.getCar_abs (2*i);
%                     %                 [distr, rel_vr] = S.getCar_abs (2*i-1);
%                     %                 if S.chkSafeDist_idxed(S.envCars(1,2),distf, rel_vf+S.envCars(1,2), 1,1) && S.chkRearSafe_idxed(S.envCars(1,2),distr, S.envCars(1,2)+rel_vr,1)
%                     if (rel_vf+S.envCars(1,2)) > maxV
%                         maxV = rel_vf+S.envCars(1,2);
%                         distY = abs(i-3)*4 +2;
%                     end
%                     %                 end
%                 end
%                 
%             end
%             if S.plot_in_loop_flag
%                 figure(1)
%                 delta = 30;
%                 plot_x = S.envCars(1,1);
%                 plot_y = distY;
%                 speed_y = -3;
%                 width = S.car_max_spd; % whatever
%                 height = 4; % whatever...
%                 xLeft = plot_x - width/2 ;
%                 yBottom = speed_y - height/2;
%                 rectangle('Position', [xLeft+S.car_min_spd, yBottom, width-S.car_min_spd, height], 'EdgeColor', 'b', 'FaceColor', [1, 0, 0, 0], 'LineWidth', 1);
%                 rectangle('Position', [xLeft+maxV-0.5, yBottom, 1, height], 'EdgeColor', 'b', 'FaceColor', [0, 1, 0, 1], 'LineWidth', 1);
%                 text(xLeft+S.car_max_spd+2,yBottom+2,num2str(maxV,3));
%                 rectangle('Position', [xLeft+S.car_min_spd, yBottom, S.envCars(1,2)-S.car_min_spd+0.1, height], 'EdgeColor', 'b', 'FaceColor', [1, 0, 0, 0.3], 'LineWidth', 1);
%                 text(S.envCars(1,2)-5,-3,num2str(S.envCars(1,2),3));
%                 rectangle('Position', [plot_x-100, plot_y-2, 300, 4], 'EdgeColor', [0, 1, 0, 0.0], 'FaceColor', [0, 1, 0, 0.1], 'LineWidth', 1);
%                 
%                 x_spd_ary = 0:0.5:10;
%                 index = 1;
%                 
%                 spd_rew_ary = 4 * exp(- 1./(x_spd_ary.^2+1)) - 2;
%                 
%                 spd_rew = 4 * exp(- 1./(observation(5).^2+1)) - 2;
%                 
%                 figure(2);clf; hold on;
%                 plot(x_spd_ary,spd_rew_ary);
%                 plot(observation(5), spd_rew,'r*');
%                 title('spd rew');
%                 figure(3);clf; hold on;
%                 
%                 y_pos_ary = 0:12;
%                 y_near_obs = observation(2)+ S.envCars(1,3);
%                 
%                 index = 1;
%                 for y_pose_temp = y_pos_ary
%                     if (y_pose_temp <= 10) && (y_pose_temp >= 2)
%                         rewY_ary(index) = 5 * exp(- 1./((y_near_obs-y_pose_temp).^2+1));
%                     elseif (y_pose_temp < 2)
%                         rewY_ary(index) = 5 * exp(- 1./((y_near_obs-2).^2+1)) -2 * abs(abs(y_pose_temp)-2) ;
%                     elseif (y_pose_temp > 10)
%                         rewY_ary(index) = 5 * exp(- 1./((y_near_obs-10).^2+1)) - 2 * abs(abs(y_pose_temp)-10) ;
%                     end
%                     index = index +1;
%                 end
%                 if (S.envCars(1,3) <= 10) && (S.envCars(1,3) >= 2)
%                     rewY = 5 * exp(- 1./((y_near_obs-S.envCars(1,3)).^2+1));
%                 elseif (S.envCars(1,3) < 2)
%                     rewY = 5 * exp(- 1./((y_near_obs-2).^2+1)) - 2 * abs(abs(S.envCars(1,3))-2) ;
%                 elseif (S.envCars(1,3) > 10)
%                     rewY = 5 * exp(- 1./((y_near_obs-10).^2+1)) - 2 * abs(abs(S.envCars(1,3))-10) ;
%                 end
%                 plot(y_pos_ary,rewY_ary);
%                 plot(S.envCars(1,3), rewY,'r*');
%                 title('y rew')
%                 spd_rew_ary = exp(-(x_spd_ary - maxV) .^ 2 ./ (2 .* 2 .* 2 .^ 2));
%                 spd_rew = exp(-(S.envCars(1,2) - maxV) .^ 2 ./ (2 .* 2 .* 2 .^ 2));
%                 if S.debug_figure
%                     figure(2);clf; hold on;
%                     plot(x_spd_ary,spd_rew_ary);
%                     plot(S.envCars(1,2), spd_rew,'r*');
%                     title('spd rew');
%                     figure(3);clf; hold on;
%                     y_pos_ary = 0:12;
%                     rewY_ary = exp(-(y_pos_ary- distY) .^ 2 ./ (2 * 2 ^2)) ;
%                     rewY = exp(-(S.envCars(1,3)- distY) .^ 2 ./ (2 * 2 ^2)) ;
%                     plot(y_pos_ary,rewY_ary);
%                     plot(S.envCars(1,3), rewY,'r*');
%                     title('y rew')
%                     
%                 end
%                 
%             end
%             rewV = exp(-(S.envCars(1,2) - maxV) ^ 2 / (2 * 2 * 2 ^ 2)) ;
%             rewY = exp(-(S.envCars(1,3) - distY) ^ 2 / (2 * 2 ^2)) ;
%             rewA = 1;
%             
%             
%             
%             if inLane && dist>midDist && S.envCars(1,2) < 31.8 && S.envCars(1,2) < maxV && ~actAccl
%                 rewA = 0;
%             end
%             if dist<midDist && actAccl
%                 rewA = 0;
%             end
%             
%             %rewA = 0; %very manual in reward function.
%             total_rew = (rewA/3+rewV/3 + rewX/3  +rewY/2)/3*4 ;
            
        end
        
        function [forward_dist, forward_spd] = getFrontCar(S,i)
            forward_dist = S.car_max_vis_dist;
            forward_spd  = S.car_max_spd;
            ego_x  = S.envCars(i,1);
            %             ego_v  = S.envCars(i,2);% get spd of ego
            %             ego_ln = S.determine_ego_lane(i);
            %[~,lane_num]=min([abs(ego_ln-0),abs(ego_ln-1),abs(ego_ln-2)])-1;% 0,1,2
            %find the closest car in front
            cars = S.envCars;
            cars(i,:) = [];
            same_lane_cars = cars(abs(cars(:,3) - S.envCars(i,3))<2,:);
            cars_infront = same_lane_cars(same_lane_cars(:,1) > ego_x,:);
            [~,idx] = min(abs(cars_infront(:,1) - ego_x));
            front_car = cars_infront(idx,:);
            if ~isempty(front_car)
                forward_dist = front_car(1)-ego_x;
                forward_spd = front_car(2);
            end
        end
        %1   2
        %3 e 4
        %5   6
        function [r_signed_dist, r_rel_spd] = getCar(S,i)
            ego_cur_lane = S.determine_ego_lane(2);
            ego_des_lane = S.A.lane_des;
            lane_info = [ego_cur_lane;ego_des_lane];
            r_signed_dist  = ones(2,1);
            r_rel_spd        = ones(2,1);
            for j = 1: 2
                cars = S.envCars(2:end,:);
                ego_x  = S.envCars(1,1);
                ego_v  = S.envCars(1,2);% get spd of ego
                ego_ln = lane_info(j);
                %[~,lane_num]=min([abs(ego_ln-0),abs(ego_ln-1),abs(ego_ln-2)])-1;% 0,1,2
                %find the closest car in front
                if ego_ln == 0 && (i == 5 || i ==6)
                    if i == 5
                        signed_dist = -S.car_max_vis_dist;
                        rel_spd     = S.car_min_spd - ego_v;
                        
                    else% i ==6
                        signed_dist = S.car_max_vis_dist;
                        rel_spd     = S.car_max_spd - ego_v;
                        
                    end
                    
                elseif ego_ln == 2 && (i == 1 || i == 2)
                    
                    if i == 1
                        signed_dist = -S.car_max_vis_dist;
                        rel_spd     = S.car_min_spd - ego_v;
                        
                    else% i ==2
                        signed_dist = S.car_max_vis_dist;
                        rel_spd     = S.car_max_spd - ego_v;
                        
                    end
                    
                else % it is a actual car that is being requested
                    if i == 1
                        lane_cars = cars(cars(:,5) == ego_ln + 1,:);
                        cars = lane_cars(lane_cars(:,1) < ego_x,:);
                    elseif i == 2
                        lane_cars = cars(cars(:,5) == ego_ln + 1,:);
                        cars = lane_cars(lane_cars(:,1) > ego_x,:);
                    elseif i == 3
                        lane_cars = cars(cars(:,5) == ego_ln,:);
                        cars = lane_cars(lane_cars(:,1) < ego_x,:);
                    elseif i == 4
                        lane_cars = cars(cars(:,5) == ego_ln,:);
                        cars = lane_cars(lane_cars(:,1) > ego_x,:);
                    elseif i == 5
                        lane_cars = cars(cars(:,5) == ego_ln-1,:);
                        cars = lane_cars(lane_cars(:,1) < ego_x,:);
                    elseif i == 6
                        lane_cars = cars(cars(:,5) == ego_ln-1,:);
                        cars = lane_cars(lane_cars(:,1) > ego_x,:);
                    end
                    [~,idx] = min(abs(cars(:,1) - ego_x));
                    wanted_car = cars(idx,:);
                    if ~isempty(wanted_car)
                        signed_dist = wanted_car(1)-ego_x;
                        rel_spd = wanted_car(2)-ego_v;
                    else
                        if i == 2 || i ==4 || i == 6
                            signed_dist = S.car_max_vis_dist;
                            rel_spd = S.car_max_spd - ego_v;
                        else
                            signed_dist = -S.car_max_vis_dist;
                            rel_spd = S.car_min_spd - ego_v;
                        end
                    end
                end
                
                
                
                
                %
                r_signed_dist(j) = signed_dist;
                r_rel_spd(j) = rel_spd;
            end
        end
        
        function setAction(S)
            for i = 1:size(S.envCars,1)
                [forward_dist, forward_spd] = S.getFrontCar(i);
                ego_v  = S.envCars(i,2);% get spd of ego
                
                if forward_dist <= S.idm_dist_1 && (ego_v - forward_spd) >= 2
                    S.envCars(i,6) = S.A_HM;
                elseif (forward_dist <= S.idm_dist_0 && (ego_v) >= S.car_min_spd) || (forward_dist <= S.idm_dist_2 && (ego_v - forward_spd) >= 2)
                    S.envCars(i,6) = S.A_BM;
                elseif forward_dist >= S.idm_dist_3  && (forward_spd - ego_v) >= 2
                    S.envCars(i,6) = S.A_AM;
                else
                    %no action, maintain
                    % S.envCars(i,4) = S.envCars(i,4) + 2;%vmax setting, so slow down
                    S.envCars(i,6) = S.A_MM;
                end
                
            end
        end
        
        function  [iniOb,LoggedSignals] =reset(S,widx,pidx)
            %             rng(4)
            rng('shuffle');

            S.icur = 1;
            S.epscur = S.epscur+1;
            
            S.envCars_history = cell(200,1);
            iniOb = zeros(5,1);
            LoggedSignals= struct;
            %         t_max = S.max_sim_time ;
            iter_max = S.max_sim_iterations ;
            S.plot_in_loop_flag = S.plot_while_running;
            S.W = S.get_world(widx) ;
            % get agent and planner
            S.A = S.get_agent(pidx) ;
            S.P = S.get_planner(pidx) ;
            
            % get agent and world ready
            S.W.reset() ;
            
            
             rng(S.epscur); % random start from here!
            
            %randomize start location for ego
            S.A.reset();
            
            
            
            % reset other cars on the road
%             S.reset_counter = S.reset_counter + 1;
%             if S.reset_counter > S.small_num_cars_iteration
%                 S.num_cars = S.min_num_cars + randi(S.max_num_cars-S.min_num_cars);
%             end
            
%             S.placeCars_blocked()
            
%             S.placeCars();
            
            S.W.setup();
            
            % get planner ready
            agent_info = S.A.get_agent_info() ;
            world_info = S.W.get_world_info(agent_info);
%             agent_info.lane = start_lane;
%             agent_info.spd = start_v;
            
            S.P.setup(agent_info,world_info) ;
            
            % check to make sure gif start is ready
            if S.save_gif
                S.start_gif = true ;
            end
            
            % initialize plot
            %             if  S.plot_while_running
            %                 S.plot(widx,pidx)
            %             end
            
            % preallocate for storing planning time spent
            S.planning_time_vec = nan(1,iter_max) ;
            
            % reset the stop counter
            S.stop_count = 0 ;
            S.stop_check_vec = false(1,iter_max) ;
            
            % reset the crash and goal checks just in case
            S.collision_check = false ;
            S.goal_check = false ;
            
            %             iniOb = agent_info.state(:,1);
            LoggedSignals.episode_result = 'initial';
            
            iniOb = S.W.get_ob(agent_info)
            
            %iniOb(7) = S.envCars(1,1)/S.A.time(end); %ego_average_vx
            %iniOb = S.get_scaled_indicator(iniOb);
            if S.plot_in_loop_flag
            S.plot()
            end
            %find something to store the states of other cars
            % start timing
            %         icur = 1 ;
            %         runtime = tic ;
            %         tstart = runtime ;
            %         tcur = toc(tstart);
            
        end
        
        
        
        %% plotting
        function plot(S,world_index,planner_index)
            S.vdisp('Plotting',3)
            
            if nargin < 3
                planner_index = 1 ;
                S.vdisp('Plotting planner 1 by default!',1) ;
                if nargin < 2
                    world_index = 1 ;
                    S.vdisp('Plotting world 1 by default!',1) ;
                end
            end
            
            fh = figure(S.figure_number) ;
            cla ; hold on ; axis equal ;
            
            % get agent, world, and planner
            A = S.get_agent(planner_index) ;
            W = S.get_world(world_index) ;
            P = S.get_planner(planner_index) ;
            
            if ~any(isinf(S.worlds{world_index}.bounds))
                axis(W.bounds+[0 0 -5 0]);
                
            end
            
            for plot_idx = S.plot_order
                switch plot_idx
                    case 'A'
                        A.plot()
                    case 'W'
                        W.plot()
                        %                         plot(W.obs_FRS(1,:),W.obs_FRS(2,:));
                    case 'P'
                        P.plot()
                        %                        plot(P.zono_exec(1,:),P.zono_exec(2,:),P.zono_color);
                        %                        plot(P.current_obstacles(1,:),P.current_obstacles(2,:));
                        
                    otherwise
                        error(['Simulator plot order is broken! Make sure ',...
                            'it is a string containing the characters ',...
                            'A (for agent), W (for world), and P (for ',...
                            'planner), in the order you want them to ',...
                            'plot (WAP is the default).'])
                end
            end
            xlim([S.A.state(1,end)-100 S.A.state(1,end)+130])%Simon's
            
            if S.save_gif
                if S.start_gif
                    % if gif saving is enabled, check that there isn't already a
                    % file in the current directory with that name
                    dir_content = dir(pwd) ;
                    file_name   = {dir_content.name} ;
                    check_name = [S.save_gif_filename,'_planner_',num2str(planner_index),'.gif'] ;
                    file_check  = any(cellfun(@(x) strcmp(check_name,x),file_name)) ;
                    if file_check
                        warning(['The current GIF filename already exists! ',...
                            'change the filename if you do not want to ',...
                            'overwrite your existing file!'])
                    end
                    if S.manually_resize_gif
                        S.vdisp(['Please resize the figure to the size you ',...
                            'want saved, or hit any key to continue.'])
                        pause
                    end
                end
                
                frame = getframe(fh) ;
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                
                filename = [S.save_gif_filename,'_planner_',num2str(planner_index),'.gif'] ;
                
                if S.start_gif
                    imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
                        'DelayTime',S.save_gif_delay_time) ;
                    S.start_gif = false ;
                else
                    imwrite(imind,cm,filename,'gif','WriteMode','append',...
                        'DelayTime',S.save_gif_delay_time) ;
                end
            end
        end
        
        function plot_at_time(S,t,planner_index,world_index)
            % method: plot_at_time(t)
            %
            % Plot the agent, world, and planner at the time t; this calls the
            % plot_at_time method for each of class.
            
            if nargin < 4
                world_index = 1 ;
            end
            
            if nargin < 3
                planner_index = 1 ;
            end
            
            if nargin < 2
                t = 0 ;
            end
            
            % get agent, world, and planner
            A = S.get_agent(planner_index) ;
            W = S.get_world(world_index) ;
            P = S.get_planner(planner_index) ;
            
            % set up hold
            if ~ishold
                hold on
                hold_check = true ;
            else
                hold_check = false ;
            end
            
            for plot_idx = S.plot_order
                switch plot_idx
                    case 'A'
                        A.plot_at_time(t)
                    case 'W'
                        W.plot_at_time(t)
                    case 'P'
                        P.plot_at_time(t)
                    otherwise
                        error(['Simulator plot order is broken! Make sure ',...
                            'it is a string containing the characters ',...
                            'A (for agent), W (for world), and P (for ',...
                            'planner), in the order you want them to ',...
                            'plot (WAP is the default).'])
                end
            end
            
            if hold_check
                hold off
            end
        end
        
        function clear_plot_data(S)
            % iterate through each agent, world, and planner, and clear the
            % plot data for each one
            for idx = 1:length(S.agents)
                clear_plot_data(S.agents{idx}) ;
            end
            
            for idx = 1:length(S.planners)
                clear_plot_data(S.planners{idx}) ;
            end
            
            for idx = 1:length(S.worlds)
                clear_plot_data(S.worlds{idx}) ;
            end
        end
        
        %% animate
        function animate(S,planner_index,world_index,save_animation_gif)
            % method: S.animate(planner_index,world_index,save_gif)
            %         S.animate(save_gif)
            %
            % Animate the agent, world, and planner for the duration given by
            % A.time. The time between animated frames is given by the simulator's
            % animation_time_discretization property. The planner and world to
            % animate are given by the planner_index and world_index inputs, so the
            % simulator plots S.worlds{world_index} and S.planners{planner_index}.
            %
            % One can also call this as S.animate(true) to save a GIF for planner 1
            % and world 1 (which is useful if there's only one planner and world
            % currently set up in the simulator).
            
            % parse input arguments
            if nargin == 2 && islogical(planner_index)
                save_animation_gif = planner_index ;
                planner_index = 1 ;
                world_index = 1 ;
            else
                if nargin < 4
                    save_animation_gif = false ;
                    start_animation_gif = false ;
                end
                
                if nargin < 3
                    world_index = 1 ;
                end
                
                if nargin < 2
                    planner_index = 1 ;
                end
            end
            
            if save_animation_gif
                start_animation_gif = true ;
                filename = S.animation_gif_setup() ;
            end
            
            % get agent
            A = S.get_agent(planner_index) ;
            
            % get world
            W = S.get_world(world_index) ;
            
            % get time
            t_vec = A.time(1):S.animation_time_discretization:A.time(end) ;
            
            % clear the active figure before animating
            if S.clear_plot_before_animating_flag
                clf ;
            end
            
            % set hold
            hold_check = ~ishold ;
            if hold_check
                hold on
            end
            
            for t_idx = t_vec
                % create plot
                S.plot_at_time(t_idx,planner_index,world_index)
                
                if S.set_plot_linewidths_flag
                    set_plot_linewidths(S.animation_linewidths) ;
                end
                
                if S.set_axes_while_animating_flag
                    axis(W.bounds)
                end
                
                % create gif
                if save_animation_gif
                    % get current figure
                    fh = get(groot,'CurrentFigure') ;
                    frame = getframe(fh) ;
                    im = frame2im(frame);
                    [imind,cm] = rgb2ind(im,256);
                    
                    if start_animation_gif
                        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
                            'DelayTime',S.animation_time_discretization) ;
                        start_animation_gif = false ;
                    else
                        imwrite(imind,cm,filename,'gif','WriteMode','append',...
                            'DelayTime',S.animation_time_discretization) ;
                    end
                else
                    pause(S.animation_time_discretization)
                end
            end
            
            if hold_check
                hold off
            end
        end
        
        function filename = animation_gif_setup(S)
            filename = S.save_gif_filename ;
            
            dir_content = dir(pwd) ;
            filenames   = {dir_content.name} ;
            file_check  = any(cellfun(@(x) strcmp(filename,x),filenames)) ;
            filename_new = filename ;
            cur_int = 1 ;
            
            while file_check
                filename_new = [filename(1:end-4),'_',num2str(cur_int),filename(end-3:end)] ;
                file_check  = any(cellfun(@(x) strcmp(filename_new,x),filenames)) ;
                cur_int = cur_int + 1 ;
            end
            
            filename = filename_new ;
        end
        
        %% utility
        function A = get_agent(S,idx)
            if nargin < 2 || S.N_agents == 1
                A = S.agents{1} ;
            else
                A = S.agents{idx} ;
            end
        end
        
        function W = get_world(S,idx)
            W = S.worlds{idx} ;
        end
        
        function P = get_planner(S,idx)
            P = S.planners{idx} ;
        end
        
        function vdisp(S,s,l,use_header)
            % Display a string 's' if the verbosity is greater than or equal to
            % the level 'l'; by default, the level is set to 0 and the default
            % threshold for displaying a message is 1 (so messages do not
            % display by default)
            if nargin < 4
                use_header = true ;
                if nargin < 3
                    l = 1 ;
                end
            end
            
            if S.verbose >= l
                if use_header
                    disp(['S: ',s])
                else
                    disp(s)
                end
            end
        end
    end
end
%         function summary = run(S,planner_indices)
%             % Method: run
%             %
%             % This function simulates the agent in the provided world as it
%             % attempts to reach the goal from its provided start position.
%
%             S.vdisp('Running simulation')
%
%             % get simulation info
%             t_max = S.max_sim_time ;
%             iter_max = S.max_sim_iterations ;
%             plot_in_loop_flag = S.plot_while_running ;
%
%             % get world and planner indices
%             world_indices = 1:length(S.worlds) ;
%
%             if nargin < 2
%                 planner_indices = 1:length(S.planners) ;
%             end
%
%             % set up summaries
%             LW = length(world_indices) ;
%             summary = cell(1,LW) ;
%
%         %% world loop
%             for widx = world_indices
%                 W = S.get_world(widx) ;
%
%                 % set up summary objects
%                 L = length(planner_indices) ;
%                 agent_name_cell = cell(1,L) ;
%                 planner_name_cell = cell(1,L) ;
%                 planner_info_cell = cell(1,L) ;
%                 agent_info_cell = cell(1,L) ;
%                 trajectory_cell = cell(1,L) ;
%                 total_real_time_cell = cell(1,L) ;
%                 total_iterations_cell = cell(1,L) ;
%                 planning_times_cell = cell(1,L) ;
%                 collision_check_cell = cell(1,L) ;
%                 goal_check_cell = cell(1,L) ;
%                 stop_check_cell = cell(1,L) ;
%                 total_simulated_time_cell = cell(1,L) ;
%                 control_input_cell = cell(1,L) ;
%                 control_input_time_cell = cell(1,L) ;
%                 t_plan_cell = cell(1,L) ;
%                 t_move_cell = cell(1,L) ;
%                 obstacles_cell = cell(1,L) ;
%
%             %% planner loop
%                 for pidx = planner_indices
%                     S.vdisp(['Planner ',num2str(pidx)])
%
%                     % get agent and planner
%                     A = S.get_agent(pidx) ;
%                     P = S.get_planner(pidx) ;
%
%                     % get agent and world ready
%                     W.reset() ;
%
%                     if isprop(A,'desired_initial_condition')
%                         A.reset([W.start; A.desired_initial_condition]) ;
%                     else
%                         A.reset(W.start) ;
%                     end
%
%
%                     % get planner ready
%                     agent_info = A.get_agent_info() ;
%                     world_info = W.get_world_info(agent_info,P) ;
%                     P.setup(agent_info,world_info) ;
%
%                     % check to make sure gif start is ready
%                     if S.save_gif
%                         S.start_gif = true ;
%                     end
%
%                     % initialize plot
%                     if plot_in_loop_flag
%                         S.plot(widx,pidx)
%                     end
%
%                     % preallocate for storing planning time spent
%                     planning_time_vec = nan(1,iter_max) ;
%
%                     % reset the stop counter
%                     S.stop_count = 0 ;
%                     stop_check_vec = false(1,iter_max) ;
%
%                     % reset the crash and goal checks just in case
%                     collision_check = false ;
%                     goal_check = false ;
%
%                     % start timing
%                     icur = 1 ;
%                     runtime = tic ;
%                     tstart = runtime ;
%                     tcur = toc(tstart);
%
%                 %% simulation loop
%                     while icur < (iter_max+1) && tcur < t_max
%                         S.vdisp('--------------------------------',3,false)
%                         S.vdisp(['ITERATION ',num2str(icur),' (t = ',num2str(A.time(end)),')'],2,false)
%
%                     %% get agent info
%                         agent_info = A.get_agent_info() ;
%
%                     %% get world info
%                         % given the current state of the agent, query the world
%                         % to get the surrounding obstacles
%                         world_info = W.get_world_info(agent_info,P) ;
%
%                     %% replan
%                         % given the current state and obstacles, query the
%                         % current planner to get a control input
%                         t_plan_spent = tic ;
%                         if S.allow_replan_errors
%                             [T_nom,U_nom,Z_nom] = P.replan(agent_info,world_info) ;
%                         else
%                             try
%                                 [T_nom,U_nom,Z_nom] = P.replan(agent_info,world_info) ;
%                             catch
%                                 S.vdisp(['Planner ',num2str(pidx),' errored while ',...
%                                          'replanning!'])
%                                 T_nom = [] ; U_nom = [] ; Z_nom = [] ;
%                             end
%                         end
%                         t_plan_spent = toc(t_plan_spent) ;
%                         planning_time_vec(icur) = t_plan_spent ;
%                         S.vdisp(['Planning time: ',num2str(t_plan_spent),' s'],4)
%
%                     %% move agent
%                         % update the agent using the current control input, so
%                         % either stop if no control was returned, or move the
%                         % agent if a valid input and time vector were returned
%                         if size(T_nom,2) < 2 || size(U_nom,2) < 2 || T_nom(end) == 0
%                             S.vdisp('Stopping!',2)
%                             A.stop(P.t_move) ;
%
%                             stop_check_vec(icur) = true ;
%
%                             % give planner a chance to recover from a stop
%                             S.stop_count = S.stop_count + 1 ;
%                             if S.stop_count > S.stop_threshold
%                                 break
%                             end
%                         else
%                             S.stop_count = 0 ;
%
%                             if ~isempty(P.t_move)
%                                 if P.t_move > T_nom(end)
%                                     S.vdisp(['The provided time vector for the ',...
%                                         'agent input is shorter than the amount of ',...
%                                         'time the agent must move at each ',...
%                                         'planning iteration. The agent will only ',...
%                                         'be moved for the duration of the ',...
%                                         'provided time vector.'],3)
%
%                                     t_move = T_nom(end) ;
%                                 else
%                                     t_move = P.t_move ;
%                                 end
%                             else
%                                 error(['Planner ',num2str(pidx),...
%                                        '''s t_move property is empty!'])
%                             end
%
%                             A.move(t_move,T_nom,U_nom,Z_nom) ;
%                         end
%
%                     %% Note (22 July 2019)
%                     % Dynamic obstacles are treated as follows:
%                     %   1) W.get_world_info should return a prediction
%                     %   2) the agent is moved according to the prediction
%                     %   3) the world moves the obstacles (according to the
%                     %      agent's movement data if needed) and then checks
%                     %      for collisions in W.collision_check
%
%                     %% crash and goal check
%                         % check if the agent is near the desired goal or if it
%                         % crashed
%                         S.vdisp('Checking if agent reached goal or crashed...',3)
%                         agent_info = A.get_agent_info() ;
%                         goal_check = W.goal_check(agent_info) ;
%                         collision_check = W.collision_check(agent_info,false) ;
%
%                         if isa(A,'multi_link_agent')
%                             S.vdisp('Checking for self-intersection.',2)
%                             collision_check = collision_check || A.self_intersection_flag ;
%                         end
%
%                         if collision_check && S.stop_sim_when_crashed
%                             S.vdisp('Crashed!',2) ;
%                             break
%                         end
%
%                         if goal_check
%                             S.vdisp('Reached goal!',2) ;
%                             break
%                         end
%
%                         % plotting and animation
%                         if plot_in_loop_flag
%                             S.plot(widx,pidx)
%                             if S.save_gif
%                                 error('Shreyas this is unfinished!')
%                             else
%                                 pause(S.plotting_pause_time) ;
%                             end
%                         end
%
%                         % pause for user if needed
%                         if S.manual_iteration
%                             user_pause = tic ;
%                             S.vdisp('Pausing for user. Press any key to continue.',2)
%                             pause
%                             user_pause = toc(user_pause) ;
%                         else
%                             user_pause = 0 ;
%                         end
%
%                         % iterate and increment time
%                         S.vdisp(['END ITERATION ',num2str(icur)],4,false)
%                         icur = icur + 1 ;
%                         tcur = toc(tstart) - user_pause ;
%                     end
%                     runtime = toc(runtime) ;
%                     S.vdisp(['Planning time spent: ',num2str(runtime)],5)
%
%                     % plot the last portion of the agent's trajectory after the
%                     % simulation ends
%                     if plot_in_loop_flag
%                         S.plot(widx,pidx)
%                     end
%
%                     S.vdisp(['Planner ',num2str(pidx), ' simulation complete!'])
%
%                 %% create summary (for the current planner)
%                     % get results at end of simulation
%                     S.vdisp('Compiling summary',3)
%                     Z = A.state ;
%                     T_nom = A.time ;
%                     U_nom = A.input ;
%                     TU = A.input_time ;
%                     agent_info = A.get_agent_info() ;
%
%                     if S.collision_check_full_traj_after_sim_flag
%                         S.vdisp('Running final collision check.',4)
%                         collision_check = W.collision_check(agent_info) ;
%                     end
%
%                     S.vdisp('Running final goal check',4)
%                     goal_check = W.goal_check(agent_info) ;
%                     agent_info_cell{pidx} = agent_info ;
%
%                     if S.save_planner_info
%                         planner_info_cell{pidx} = S.planners{pidx}.info ;
%                     else
%                         planner_info_cell{pidx} = 'no info saved' ;
%                     end
%
%                     % fill in the results for the current planner
%                     agent_name_cell{pidx} = A.name ;
%                     planner_name_cell{pidx} = P.name ;
%                     trajectory_cell{pidx} = Z ;
%                     total_simulated_time_cell{pidx} = T_nom ;
%                     control_input_cell{pidx} = U_nom ;
%                     control_input_time_cell{pidx} = TU ;
%                     total_real_time_cell{pidx} = runtime ;
%                     total_iterations_cell{pidx} = icur ;
%                     planning_times_cell{pidx} = planning_time_vec ;
%                     collision_check_cell{pidx} = collision_check ;
%                     goal_check_cell{pidx} = goal_check ;
%                     stop_check_cell{pidx} = stop_check_vec ;
%                     t_plan_cell{pidx} = P.t_plan ;
%                     t_move_cell{pidx} = P.t_move ;
%                     obstacles_cell{pidx} = W.obstacles;
%
%                     if goal_check
%                         S.vdisp('In final check, agent reached goal!')
%                     end
%
%                     if collision_check
%                         S.vdisp('In final check, agent crashed!')
%                     end
%                 end
%                 S.vdisp(['World ',num2str(widx),' complete! Generating summary.'])
%                 summary{widx} = struct('agent_name',agent_name_cell,...
%                                  'planner_name',planner_name_cell,...
%                                  'trajectory',trajectory_cell,...
%                                  'total_real_time',total_real_time_cell,...
%                                  'total_iterations',total_iterations_cell,...
%                                  'planning_time',planning_times_cell,...
%                                  'collision_check',collision_check_cell,...
%                                  'goal_check',goal_check_cell,...
%                                  'stop_check',stop_check_cell,...
%                                  'total_simulated_time',total_simulated_time_cell,...
%                                  'control_input',control_input_cell,...
%                                  'control_input_time',control_input_time_cell,...
%                                  'agent_info',agent_info_cell,...
%                                  'planner_info',planner_info_cell,...
%                                  'obstacles',obstacles_cell,...
%                                  'planner_indices',planner_indices,...
%                                  'N_obstacles',W.N_obstacles,...
%                                  't_plan',t_plan_cell,...
%                                  't_move',t_move_cell,...
%                                  't_max',t_max,...
%                                  'iter_max',iter_max,...
%                                  'start',W.start,...
%                                  'goal',W.goal,...
%                                  'bounds',W.bounds,...
%                                  'notes','') ;
%             end
%
%             % clean up summary if only one world was run
%             if LW == 1
%                 summary = summary{1} ;
%             end
%
%             S.simulation_summary = summary ;
%
%             if nargout < 1
%                 clear summary ;
%                 S.vdisp('Simulation summary stored in simulator_obj.simulation_summary.',1)
%             end
%
%             S.vdisp('Simulation complete!')
%         end
