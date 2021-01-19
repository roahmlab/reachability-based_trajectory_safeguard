classdef quadrotor_agentHelper < agentHelper
    %% properties
    properties
        % FRS info
        FRS
        FRS_time_step
        FRS_time
        FRS_N_steps
        
        % tracking error
        tracking_error_table
        tracking_error_type = 'table' ; % 'none', 'constant', or 'table'
        tracking_error_zono = [] ;
        tracking_error_constant_value = 0.1 ;
        
        % planned trajectory timing
        t_peak = 1 ; % desired time to peak speed
        t_total = 3 ; % total duration of trajectory
        t_sample = 0.01 ; % time discretization
        t_extra = 2 ; % extra time added to end of trajectory so the
        % low-level controller has a reference after the
        % quadrotor comes to a stop
        
        % property to save the current planned trajectory
        current_parameter
        
        proposed_parameter
        % how to get initial condition for previous trajectory
        use_agent_for_initial_condition_flag = true ;
        
        % speed and acceleration info
        v_max = 2 ;
        a_max = 3 ;
        a_est = zeros(3,1) ;
        
        % properties for sampling-based method of adjusting params
        N_sample = 15; % should be an odd number to ensure symmetry
        
        % plotting
        plot_local_view = false ;
        plot_zonotope_reach_set_flag = true ;
        plot_zono_style = 'tube' ; % 'face' or 'tube'
        plot_zono_face_color = 1/256*[199,234,229] ;
        plot_zono_face_opacity = 0.0 ;
        plot_zono_edge_color = 1/256*[199,234,229] ;
        plot_zono_edge_opacity = 0.5 ;
        plot_zono_edge_width = 2 ;
        plot_zono_lighting = 'flat' ;
        plot_zono_skip_idxs = 5 ; % # of time steps to skip when plotting zonos
        plot_desired_trajectory_line_width = 3 ;
        
        k_adjusted = [0;0;0];
        ref_Z= [];
        proposed_ref_Z= [];
        t_real_start = [];
        t_proposed_start = [];
        plot_flag = 0;
        FRS_plotting_param=[];
        
        S = [];
        HLP = [];
        
        fminconopt = [];
        add_noise_to_waypoint = false;
        N_stop = 0;
        N_stop_threshold = 3;
        delta_v_max = 2;
        
%         t_plan_vec = [];
        sphere_density = 10;
    end
    
    %% methods
    methods
        %% agentHelper required methods start here
        %% constructor
        function AH = quadrotor_agentHelper(A,FRS,tracking_error_table,HLP,varargin)
            % AH = quadrotor_agentHelper(A,FRS,tracking_error_table)
            %
            % Create an agentHelper for safe RL with the quadrotor
            
            % create superclass instance
            AH@agentHelper() ;
            
            % parse varargin
            AH = parse_args(AH,'name','quadrotor_agentHelper',varargin{:}) ;
            
            % get agent
            AH.A = A ;
            
            % get FRS zonotopes (the Rcont property of the zonotope FRS
            % that is loaded from a .mat file)
            AH.vdisp('Loading FRS',6)
            AH.FRS = FRS.Rcont ;
            
            % get the tracking error table
            if strcmp(AH.tracking_error_type,'table')
                AH.vdisp('Setting up tracking error table',6)
                if isstruct(tracking_error_table)
                    AH.tracking_error_table = tracking_error_table.tracking_error_table ;
                elseif istable(tracking_error_table)
                    AH.tracking_error_table = tracking_error_table ;
                else
                    error('Please pass the tracking error table in correctly!')
                end
            end
            
            % create tracking error zonotope for the constant case
            tracking_c = zeros(size(AH.FRS{1}{1}.Z, 1), 1);
            tracking_G = zeros(size(AH.FRS{1}{1}.Z, 1), 3);
            tracking_G(1, 1) = AH.tracking_error_constant_value ;
            tracking_G(6, 2) = AH.tracking_error_constant_value ;
            tracking_G(11, 3) = AH.tracking_error_constant_value ;
            AH.tracking_error_zono = zonotope([tracking_c, tracking_G]);
            
            % try extracting timing info from FRS
            try
                % extract trajectory timing
                AH.t_peak = FRS.t_peak ;
                AH.t_total = FRS.t_total ;
                AH.FRS_time_step = FRS.options1.timeStep ;
                AH.FRS_time = 0:AH.FRS_time_step:AH.t_total ;
                AH.FRS_N_steps = length(AH.FRS) ;
                
                % set agentHelper timing
                if isempty(AH.t_move) || (AH.t_move == 0)
                    AH.vdisp('Setting agentHelper''s t_move property automatically',3)
                    AH.t_move = AH.t_peak / 2 ;
                end
                AH.t_failsafe_move = AH.t_total - AH.t_move ;
                
                AH.vdisp('Timing info extracted from FRS',2)
            catch
                AH.vdisp('Timing info not included in FRS file! Using defaults.',2)
            end
            
            if ~isempty(HLP)
                AH.HLP = HLP;
            end
            
            o = optimoptions('fmincon') ;
            o.OptimalityTolerance = 1e-3 ;
            o.MaxIterations = 100000 ;
            o.MaxFunctionEvaluations = 10000 ;
            o.SpecifyConstraintGradient = true ;
            o.SpecifyObjectiveGradient = isa(AH.HLP,'highway_HLP');
            o.Display = 'off';
            o.CheckGradients = false;
    %             end
            AH.fminconopt = o;
            % set up plot_data
            AH.plot_data.FRS = [] ;
            AH.plot_data.trajectory = [] ;
        end
        
        %% reset
        function reset(AH,flags,seed)
            % quadrotor_agentHelper.reset(flags,world_start)
            if exist('seed','var')
                rng(seed)
            end
            % Reset the quadrotor agent to either the provided world start
            % or to its default, and reset the agentHelper's flags struct.
            AH.k_adjusted = [0;0;0];
            AH.vdisp('Resetting agent',4)
            %             if nargin < 3
            AH.A.reset([1;0;0]) ;
            %             else
            %                 AH.A.reset(world_start) ;
            %             end
            
            AH.vdisp('Setting flags',4)
            AH.flags = flags ;
            AH.ref_Z = [];
            AH.proposed_ref_Z = [];
            if ~isempty(AH.HLP)
                AH.HLP.waypoints = [] ;
            end
            if ~isempty(AH.S)
                if AH.S.discrete_flag
                    AH.N_sample = 5;
                    AH.sphere_density= 3;
                end
            end
        end
        
        %% gen_param: R mode generate a parameter based on a cost function
        function [k, no_replace_action] = gen_param(AH, world_info)
            replan_start_tic = tic;
            % do all the setup stuff required for replanning
            agent_info = AH.get_agent_info() ;
            [x_0,v_0,a_0] = AH.get_initial_condition(agent_info) ;
            
            % get zonotopes from world and shift them to the agent's
            % current position
            O_zono = AH.process_obstacles(world_info.obstacles,x_0) ;
            
            % get desired waypoint IN GLOBAL COORDINATES
            lookahead_distance = AH.t_peak*AH.v_max ;
            if ~isempty(AH.HLP)
                x_des = AH.HLP.get_waypoint(agent_info,world_info,lookahead_distance) ;
                if AH.add_noise_to_waypoint
                mn = 0 ;
                sd = 0.5*(7/(1+norm(v_0))^4) ;
                x_des = x_des + [rand_range(-sd,sd,mn,sd,2,1); 0] ;
                end
                x_des = x_des - x_0 ;
            else
                x_des = [];% no hlp, optimize reward directly
            end
                try
                [A_con,b_con] = AH.generate_constraints(x_0,v_0,a_0,O_zono) ;

%                 [k,~] = AH.trajopt_fmincon(A_con,b_con,v_0,a_0,x_des,replan_start_tic) ;
                [k,exitflag] =AH.trajopt_sample(A_con,b_con,v_0,a_0,x_0,x_des,replan_start_tic) ;
                    if exitflag ~= 1 || isempty(k)
                        no_replace_action =1;
                    else
                        no_replace_action =0;
                    end
                catch
    %                     P.vdisp('Unable to find new trajectory!',3)
                    k = nan(3,1) ;
                    no_replace_action =1;

                end
            
             
%               AH.t_plan_vec = [AH.t_plan_vec toc(replan_start_tic)];
            
                % FOR DEBUGGING:
                % plot3(x_des(1),x_des(2),x_des(3),'ko')
%             end
            
            % move x_des to local coordinates
            
%             no_replace_action = 0;
            
        % check if the previous plan 
%                 Ah.vdisp('Replanning',3)
            
            
            
            
            % create dummy control input
%             N_t_des = size(T_des,2) ;
%             U_des = zeros(4,N_t_des) ;
%             
%             % update acceleration estimate
%             a_est_new = match_trajectories(P.t_plan,T_des,Z_des(7:9,:)) ;
%             
%             % update current plan and obstacles
%             P.current_plan.T_des = T_des ;
%             P.current_plan.Z_des = Z_des ;
%             P.current_plan.a_est = Z_des(7:9,:) ;
%             
%             % update planner info
%             P.info.obstacles = [P.info.obstacles, {P.current_obstacles}] ;
%             P.info.t_start_plan = [P.info.t_start_plan,agent_info.time(end)] ;
%             P.info.x_0 = [P.info.x_0,x_0] ;
%             P.info.v_0 = [P.info.v_0,v_0] ;
%             P.info.a_0 = [P.info.a_0,a_0] ;
%             P.info.v_peak = [P.info.v_peak,v_peak] ;
%             P.info.a_est = [P.info.a_est, a_est_new] ;
%             P.info.T_des = [P.info.T_des, {T_des}] ;
%             P.info.U_des = [P.info.U_des, {U_des}] ;
%             P.info.Z_des = [P.info.Z_des, {Z_des}] ;
        end
       function [v_peak,exitflag] = trajopt_sample(AH,A_con,b_con,v_0,a_0,x_0,x_des,replan_start_tic)
            % create sphere at v_0
%             sample_start = tic;
            time_limit = AH.t_peak;
            N_s = AH.N_sample ;
            % check if the quadrotor has stopped too many times and reduce
            % the number of samples if so (this can help find solutions
            % faster when stuck around too many obstacles to plan quickly)
            if norm(v_0) < 1e-3
                AH.N_stop = AH.N_stop + 1 ;
                if AH.N_stop > AH.N_stop_threshold
                    N_s = 3 ;
                end
            else
                AH.N_stop = 0 ;
            end

            S_v = make_v_peak_sphere(AH.delta_v_max,N_s,v_0,AH.sphere_density) ;
            
            error_if_out_of_time(replan_start_tic,time_limit)
            
            % remove v > v_max
            S_v_log = vecnorm(S_v) <= AH.v_max ;
            S_v = S_v(:,S_v_log) ;
            
            error_if_out_of_time(replan_start_tic,time_limit)
            
            % figure(3) ; cla ; hold on ; axis equal

            % evaluate constraints
            if ~isempty(A_con)
                C_eval = A_con*S_v + b_con ;
                C_rshp = reshape(C_eval,6,[]) ;
                C_min = min(C_rshp,[],1) + 1e-6;
                C_rshp_2 = reshape(C_min,size(C_eval,1)/6,[]) ;
                C_max = max(C_rshp_2,[],1) ;
                C_log = C_max < 0 ;
                S_v = S_v(:,C_log) ;
            end
%             display('constraints eval took')
%             toc(sample_start)% time took for constraint eval
%             sample_start = tic;
            error_if_out_of_time(replan_start_tic,time_limit)
            
            % plot3(S_v(1,:),S_v(2,:),S_v(3,:),'b.') ;
            
            % evaluate cost (note that t_peak is hard-coded in here)
            if ~isempty(x_des)
                p_peak = pos_quadrotor_peak(v_0,a_0,S_v);
                J_vals = vecnorm(p_peak - x_des) ;
            else
                p_peak = pos_quadrotor_peak(v_0,a_0,S_v) + x_0;
                aginfo = AH.get_agent_info;
                J_vals = ones(size(p_peak,2),1)*inf;
                for future_pos_idx = 1:size(p_peak,2)
                    aginfo.velocity = S_v(:,future_pos_idx);
                    aginfo.position = p_peak(:,future_pos_idx);
                    ob = AH.S.W.get_ob(aginfo);
                    J_vals(future_pos_idx) = -AH.S.W.getRew(aginfo,ob);
%                     [~, min_idx] = min(J_vals);
                    if toc(replan_start_tic) > 0.95*time_limit
                        break;
                    end
                end
            end
            [~,min_idx] = min(J_vals) ;
%             display('cost eval took')
%             toc(sample_start) %time took for cost eval
            v_peak = S_v(:,min_idx);
            exitflag = 1;
            if norm(v_peak) < 1e-3
                AH.N_stop = AH.N_stop + 1 ;
                if AH.N_stop > AH.N_stop_threshold
                   exitflag = -1;
                end
            else
                AH.N_stop = 0 ;
            end
            % FOR DEBUGGING:
            % plot3(v_peak(1),v_peak(2),v_peak(3),'ko','Markersize',12)
            error_if_out_of_time(replan_start_tic,time_limit)
       end 
        
       function [v_peak,exitflag] = trajopt_fmincon(AH,A_con,b_con,v_0,a_0,x_des,replan_start_tic)
%             P.vdisp('Setting up cost and constraint functions',7)
            % make cost and nonlcon funcs
            tlimit = AH.t_peak;
%             tlimit = 5;
            cost = @(v_peak) eval_cost(v_peak, v_0, a_0, x_des, replan_start_tic, tlimit) ;
            % cons = @(v_peak) eval_zono_cons(v_peak, A_con, b_con) ;
            cons = @(v_peak) eval_cons(v_peak, v_0, AH.v_max, AH.a_max,...
                AH.t_peak, A_con, b_con,...
                replan_start_tic, tlimit) ;
            
            % call fmincon
%             P.vdisp('Calling fmincon',7)
            epsilon = 1e-3 ;
            lb = (-5 + epsilon)*ones(3,1) ;
            ub = (5 - epsilon)*ones(3,1) ;
            
            % make initial guess towards x_des
            if vecnorm(x_des) > 0
                initial_guess = 0.25.*make_unit_length(x_des) ;
            else
                initial_guess = zeros(3,1) ;
            end
            
            
            [v_peak,~,exitflag,~] = fmincon(cost, initial_guess, [], [],...
                [], [], lb, ub, cons, AH.fminconopt) ;
        end
        %% adjust
        function [k, replace_distance, replaced_flag] = adjust(AH,k_user,world_info)
            % [k,replace_distance,replaced_flag] = AH.adjust(k_user,world_info)
            %
            % Given a choice of trajectory parameter k_user, and world
            % information (i.e., obstacles), determine if k_user is unsafe
            % and adjust it to a safe choice of k instead. Return the
            % adjusted k, the distance from the original k, and a flag that
            % states if the adjustment was successful or not (adjustment
            % may be unsuccessful when there are no safe trajectories
            % available)
            
            % save Z user for plotting as well
            real_ref_flag = 0;
            AH.gen_ref(k_user, real_ref_flag);
            
            
            % display the user's desired velocity
            AH.vdisp(['Desired velocity (x,y,z): (',...
                num2str(k_user(1)),',',...
                num2str(k_user(2)),',',...
                num2str(k_user(3)),')'],5)
            
            % get agent info and initial condition
            agent_info = AH.get_agent_info() ;
            [x_0,v_0,a_0] = AH.get_initial_condition(agent_info) ;
            
            % initialize replace distance to 0 optimistically; this will be
            % nonzero if k_user is chosen to be unsafe
            replace_distance = 0 ;
            replaced_flag = 0 ;
            
            % check that the agent initial condition is sliceable
            if (vecnorm(v_0) >= AH.v_max)
                v_0 = make_unit_length(v_0).*AH.v_max ;
                warning('Fixed agent initial velocity!')
            end
            
            % get obstacles from world and shift them to the agent's
            % current position
            AH.vdisp('Getting obstacles',4) ;
            O_zono = AH.process_obstacles(world_info.obstacles,x_0) ;
            
            % convert obstacles to constraints on k
            [A_con,b_con] = AH.generate_constraints(x_0,v_0,a_0,O_zono) ;
            
            % check if k_user is safe WRT obstacles, max speed, and
            % acceleration constraint
            k_user_check = (vecnorm(k_user) <= AH.v_max) && ...
                ((vecnorm(k_user - v_0)/AH.t_peak) <= AH.a_max) ;
            if ~isempty(A_con)
                C_eval = A_con*k_user + b_con ;
                C_rshp = reshape(C_eval,6,[]) ;
                C_min = min(C_rshp,[],1) + 1e-6;
                C_rshp_2 = reshape(C_min,size(C_eval,1)/6,[]) ;
                C_max = max(C_rshp_2,[],1) ;
                C_log = C_max < 0 ;
                k_user_check = k_user_check && (C_log) ;
            end
            
            % adjust the k_user input if it is unsafe
            if k_user_check
                AH.vdisp('Parameter does not need to be adjusted',6) ;
                k = k_user ;
            else
                AH.vdisp('Unsafe parameter entered! Trying to find a safe one',8)
                replaced_flag = 1 ;
                
                % generate a sphere of samples of peak velocity centered at
                % the drone's current velocity
                k_sphere = make_v_peak_sphere(AH.a_max.*AH.t_peak,AH.N_sample,v_0) ;
                
                % remove all choices of v > v_max
                k_sphere_log = vecnorm(k_sphere) <= AH.v_max ;
                k_sphere = k_sphere(:,k_sphere_log) ;
                
                % evaluate constraints
                if ~isempty(A_con) && ~isempty(k_sphere)
                    C_eval = A_con*k_sphere + b_con ;
                    C_rshp = reshape(C_eval,6,[]) ;
                    C_min = min(C_rshp,[],1) + 1e-6;
                    C_rshp_2 = reshape(C_min,size(C_eval,1)/6,[]) ;
                    C_max = max(C_rshp_2,[],1) ;
                    C_log = C_max < 0 ;
                    k_sphere = k_sphere(:,C_log) ;
                end
                
                % check if there are any feasible choices of k left
                if isempty(k_sphere)
                    AH.vdisp('Parameter adjustment failed!',6)
                    k = [0;0;0] ;
                    replaced_flag = 2 ;
                    replace_distance = [] ;
                else
                    AH.vdisp('Parameter adjustment successful!',6)
                    
                    % pick the closest choice to the user input
                    delta_v = vecnorm(k_sphere - repmat(k_user,1,size(k_sphere,2))) ;
                    [replace_distance,best_idx] = min(delta_v) ;
                    k = k_sphere(:,best_idx) ;
                    
                    % display the user's desired velocity
                    AH.vdisp(['Adjusted velocity (x,y,z): (',...
                        num2str(k(1)),',',...
                        num2str(k(2)),',',...
                        num2str(k(3)),')'],5)
                    if norm(k) < 1e-2
                        AH.N_stop = AH.N_stop + 1 ;
                        if AH.N_stop > AH.N_stop_threshold
                            replaced_flag = 2; %stuck if adjusted is 0 multiple times
                        end
                    else
                        AH.N_stop = 0 ;
                    end

                end
            end
            AH.k_adjusted = k;
        end
        
        %% generate reference trajectory
        function [T,U,Z] = gen_ref(AH,v_peak, real_reference_flag)
            % get agent info and initial condition
            agent_info = AH.get_agent_info() ;
            [x_0,v_0,a_0] = AH.get_initial_condition(agent_info) ;
            
            % compute desired trajectory with output
            [T,Z,~] = generate_spline_peak_speed(v_0,a_0,v_peak,...
                AH.t_move,AH.t_peak,AH.t_total,AH.t_sample,...
                AH.t_extra) ;
            
            % translate reference trajectory to where quadrotor is
            Z(1:3,:) = Z(1:3,:) + x_0 ;
            
            % create dummy control input
            N_t = length(T) ;
            U = zeros(4,N_t) ;
            
            % updated the saved reference
            if ~exist('real_reference_flag','var')
                real_reference_flag = 1;
            end
            if real_reference_flag
%                 AH.ref_Z=[AH.ref_Z;x_cur+Z(1,:);y_cur+Z(2,:)];% for plotting
                AH.t_real_start = [AH.t_real_start;AH.A.time(end)];
                AH.current_parameter = [v_0 ; a_0 ; v_peak] ;
                AH.ref_Z=[AH.ref_Z;Z(1:3,:)];% for plotting
            else
                AH.t_proposed_start = [AH.t_proposed_start;AH.A.time(end)];
                AH.proposed_parameter = [v_0 ; a_0 ; v_peak] ;
                AH.proposed_ref_Z=[AH.proposed_ref_Z;Z(1:3,:)];% for plotting
            end
        end
        
        %% convert action to parameter
        function k = convert_action_to_parameter(AH,action,~)
            % k = convert_action_to_parameter(AH,action)
            %
            % The action is an output of the RL agent in the range [-1,1]^n
            % where n is the appropriate dimension; in this case, n = 3
            % because the action is a choice of v_peak for the drone.
            % action stands for delta vx, and parametera stands for peak v
            % to convert add
            para = action + AH.k_adjusted;
            if vecnorm(para) > AH.v_max
                k = AH.v_max.*make_unit_length(para) ;
            else
                k = para;
            end
            if AH.flags.safety_layer == 'N'
                AH.k_adjusted = k; % consider putting that in AH parent
            end
        end
        
        %% plotting
        function plot(AH)
            % quadrotor_agentHelper.plot() ;
            %
            % Plot the current trajectory and corresponding FRS zonotopes.
            
            % set hold on if it isn't already on
            hold_check = hold_switch() ;
            
            % plot the current trajectory and proposed
            proposed_and_adjusted_flag = 1;
            AH.plot_current_trajectory(proposed_and_adjusted_flag);
            
            % plot the reachable set for the current trajectory
            if AH.plot_zonotope_reach_set_flag
                Z = AH.ref_Z ;
                k = AH.current_parameter ;
                
                if ~isempty(Z) && ~isempty(k)
                    x_0 = Z(1:3,1) ;
                    v_0 = k(1:3) ;
                    a_0 = k(4:6) ;
                    v_peak = k(7:9) ;
                    
                    AH.plot_zonotope_reach_set(v_peak,x_0,v_0,a_0) ;
                end
            end
            if ~isempty(AH.HLP)
                AH.HLP.plot();
            end
            % set hold off if hold wasn't on originally
            hold_switch(hold_check) ;
        end
        
        function plot_adjust(AH)
            AH.plot() ;
        end
        
        %% helper methods start here
        %% get initial condition
        function [x_0,v_0,a_0] = get_initial_condition(AH,agent_info)
            if AH.use_agent_for_initial_condition_flag || isempty(AH.T)
                % get the initial condition for the next traj
                x_0 = agent_info.position(:,end) ;
                v_0 = agent_info.velocity(:,end) ;
            else
                z_0 = match_trajectories(AH.t_move,AH.T,AH.Z) ;
                x_0 = z_0(1:3) ;
                v_0 = z_0(4:6) ;
            end
            
            % get the acceleration
            if ~isempty(AH.T)
                a_0 = match_trajectories(AH.t_move,AH.T,AH.Z(7:9,:)) ;
            else
                a_0 = AH.a_est(:,end) ;
            end
        end
        
        %% process obstacles
        function O_zono = process_obstacles(AH,O,x_0)
            % get obstacle number
            N_obstacles = length(O) ;
            AH.vdisp(['N obstacles: ',num2str(length(O))],5) ;
            
            % fill in a cell array that is the length of the FRS in time
            % steps, where each cell contains all a cell array of all the
            % zonotope obstacles at that time
            O_zono = cell(1,AH.FRS_N_steps) ;
            
            for FRS_idx = 1:AH.FRS_N_steps
                % create cell array to put into full cell array that will
                % contain each obstacle zonotope for the current time step
                O_zono_idx = cell(1,N_obstacles) ;
                
                for o_idx = 1:N_obstacles
                    % get the current obstacle at the current time step
                    o = O{o_idx} ;
                    
                    if isa(o,'zonotope_dynamic_obstacle')
                        error('dynamic obstacles not yet supported!')
                    else
                        Z = o.zono - x_0 ;
                    end
                    
                    % fill in the obstacle cell array
                    O_zono_idx{o_idx} = Z ;
                end
                
                % fill the FRS obstacles in
                O_zono{FRS_idx} = O_zono_idx ;
            end
        end
        
        %% convert obstacles to constraints on k
        function [A_con, b_con] = generate_constraints(AH,~,v_0,a_0,O_zono)
            % generate constraint matrices for trajopt
            switch AH.tracking_error_type
                case 'none'
                    [A_con, b_con] = generate_quadrotor_trajopt_constraints(v_0, a_0,...
                        AH.FRS, O_zono) ;
                case 'constant'
                    [A_con, b_con] = generate_quadrotor_trajopt_constraints(v_0, a_0,...
                        AH.FRS, O_zono, AH.tracking_error_zono) ;
                case 'table'
                    [A_con, b_con] = generate_quadrotor_trajopt_constraints(v_0, a_0,...
                        AH.FRS, O_zono, AH.tracking_error_table) ;
                otherwise
                    error('Tracking error type must be none, constant, or table.')
            end
        end
        
        %% plotting helper methods
        %% plot FRS
        function plot_zonotope_reach_set(AH,v_pk,varargin)
            % AH.plot_zonotope_reach_set(v_pk)
            % AH.plot_zonotope_reach_set(v_pk,x_0,v_0,a_0)
            %
            % Plot the zonotope reach set for the provided v_peak value and
            % (optionally) the provided initial conditions
            
            if nargin < 3
                [x_0,v_0,a_0] = get_initial_condition(AH,AH.A.get_agent_info()) ;
            else
                x_0 = varargin{1} ;
                v_0 = varargin{2} ;
                a_0 = varargin{3} ;
            end
            
            [F,V] = AH.make_zonotope_reach_set_patch_data(v_pk,x_0,v_0,a_0) ;
            
            if check_if_plot_is_available(AH,'FRS')
                AH.plot_data.FRS.Faces = F ;
                AH.plot_data.FRS.Vertices = V ;
            else
                patch_data = patch('Faces',F,...
                    'Vertices',V,...
                    'FaceColor',AH.plot_zono_face_color,...
                    'FaceAlpha',AH.plot_zono_face_opacity,...
                    'EdgeColor',AH.plot_zono_edge_color,...
                    'EdgeAlpha',AH.plot_zono_edge_opacity,...
                    'LineWidth',AH.plot_zono_edge_width) ;
                AH.plot_data.FRS = patch_data ;
            end
        end
        
        %% make FRS plotting object
        function [F_out,V_out] = make_zonotope_reach_set_patch_data(AH,v_pk,x_0,v_0,a_0)
            % get the FRS for plotting
            FRS_plot = AH.FRS ;
            
            % add tracking error
            switch AH.tracking_error_type
                case 'constant'
                    FRS_plot = add_tracking_error_to_FRS(FRS_plot, AH.tracking_error_zono);
                case 'table'
                    FRS_plot = add_tracking_error_to_FRS(FRS_plot, AH.tracking_error_table, v_0);
            end
            
            % slice the FRS by k = (k_pk,k_v,k_a)
            FRS_plot = slice_FRS(FRS_plot, [2; 7; 12; 3; 8; 13; 4; 9; 14], [v_0; a_0; v_pk]);
            
            % generate patch info for the FRS
            position_dimensions = [1; 6; 11];
            F_out = [] ;
            V_out = [] ;
            cnt = 0;
            for i = 1:AH.plot_zono_skip_idxs:length(FRS_plot)
                Z = FRS_plot{i}{1}.Z;
                c = Z(position_dimensions, 1) + x_0;
                G = Z(position_dimensions, 2:end);
                lwh = 2*sum(abs(G),2);
                [F,V] = make_cuboid_for_patch(lwh(1),lwh(2),lwh(3),c) ;
                
                switch AH.plot_zono_style
                    case 'face'
                        % get the front face and corresponding vertices
                        F = F(4,:) ;
                        V = V(F(:),:) ;
                        F_out = [F_out ; (4*cnt) + (1:4)] ;
                    case 'tube'
                        % keep all faces
                        F_out = [F_out ; (8*cnt) + F] ;
                end
                V_out = [V_out ; V] ;
                cnt = cnt + 1 ;
            end
        end
        
        %% plot trajectory
        function plot_current_trajectory(AH,proposed_and_adjusted)
            % AH.plot_current_trajectory(t)
            %
            % Plot the current trajectory (which is in the
            % AH.current_trajectory property)
            
%             if ~isempty(AH.proposed_ref_Z)
%                 plot(AH.proposed_ref_Z(end-1,:),AH.proposed_ref_Z(end,:),'k-','LineWidth',3);
%                 plot(AH.proposed_ref_Z(end-1,:),AH.proposed_ref_Z(end,:),'Color','y','LineWidth',3,'LineStyle','--');
%                 
%                 %                 xlim([AH.ref_Z(1,1)-20,AH.ref_Z(1,1)+30]);
%             end
%             if ~isempty(AH.ref_Z)
%                 plot(AH.ref_Z(end-1,:),AH.ref_Z(end,:),'Color',[0 0 0],'LineStyle','-','LineWidth',3);
%                 %                 plot(AH.ref_Z(end-1,:),AH.ref_Z(end,:),'g--','LineWidth',2);
%                 
%                 %                 xlim([AH.ref_Z(1,1)-20,AH.ref_Z(1,1)+30]);
%             end
            
            if proposed_and_adjusted
                if ~isempty(AH.proposed_ref_Z)
                    Z = AH.proposed_ref_Z(end-2:end,:) ;
                    
                    % plot current plan
                    if check_if_plot_is_available(AH,'proposed_ref_Z')
                        AH.plot_data.proposed_ref_Z.XData = Z(1,:) ;
                        AH.plot_data.proposed_ref_Z.YData = Z(2,:) ;
                        AH.plot_data.proposed_ref_Z.ZData = Z(3,:) ;
                        AH.plot_data.proposed_trajectory_overlay.XData = Z(1,:) ;
                        AH.plot_data.proposed_trajectory_overlay.YData = Z(2,:) ;
                        AH.plot_data.proposed_trajectory_overlay.ZData = Z(3,:) ;
                    else
                        %                     trajectory_data = plot3(Z(1,:),Z(2,:),Z(3,:),':k',...
                        %                         'LineWidth',AH.plot_desired_trajectory_line_width) ;
                        trajectory_data = plot3(Z(1,:),Z(2,:),Z(3,:),'k-','LineWidth',AH.plot_desired_trajectory_line_width);
                        trajectory_data2 = plot3(Z(1,:),Z(2,:),Z(3,:),'Color','y','LineWidth',AH.plot_desired_trajectory_line_width,'LineStyle','--');
                        AH.plot_data.proposed_ref_Z = trajectory_data ;
                        AH.plot_data.proposed_trajectory_overlay = trajectory_data2 ;
                    end
                end
            end
            if ~isempty(AH.ref_Z)
                Z = AH.ref_Z(end-2:end,:)  ;
                
                % plot current plan
                if check_if_plot_is_available(AH,'trajectory')
                    AH.plot_data.ref_Z.XData = Z(1,:) ;
                    AH.plot_data.ref_Z.YData = Z(2,:) ;
                    AH.plot_data.ref_Z.ZData = Z(3,:) ;
                else
                    trajectory_data = plot3(Z(1,:),Z(2,:),Z(3,:),'k-',...
                        'LineWidth',AH.plot_desired_trajectory_line_width) ;
                    AH.plot_data.ref_Z = trajectory_data ;
                end
            end
            
            
            
        end
    end
end



