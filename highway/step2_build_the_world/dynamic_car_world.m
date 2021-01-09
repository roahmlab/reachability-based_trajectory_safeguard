classdef dynamic_car_world < world
    properties
        buffer = 1 ;
        obstacles_unseen
        obstacles_seen
        obstacle_size_bounds = [0.2 0.6];
        obstacle_rotation_bounds = [0,0] ;
        bounds_as_polyline
        
        % plotting
        obstacle_seen_color = [1 0 0] ;
        obstacle_unseen_color = [1 0.6 0.6] ;
        obstacle_history
        obstacle_alpha = 0.5;
        
        obstacle_size = [4.8, 2] ;
        %         name = 'dyn_obs_world'
        envCars
        num_cars = 10;
        
        o %this is obstacle centered at 0
        o_no_frs
        v_array = [22 24 26 28 30 32];
        obs_FRS
        obstacles_old % for collision check
        
        SIM_MAX_DISTANCE = 500;
        car_safe_dist = 26;%40 % car 6 m long
        car_max_vis_dist =120;
        
        car_max_spd = 5;
        car_min_spd = 1;
        
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
        start_line
        
    end
    
    methods
        %% constructor
        function W = dynamic_car_world(varargin)
            % set default properties
            start_line = varargin{2}(1);
            bounds = [start_line -start_line 0 12] ;
            W.start_line = -start_line;
            N_obstacles = 0 ;
            
            % parse other args
            W = parse_args(W,'bounds',bounds,'N_obstacles',N_obstacles,varargin{:}) ;
            
            % call setup to generate obstacles
            W.setup() ;
            
            % check for mapping toolbox
            try
                polyxpoly([],[],[],[]) ;
            catch
                error(['Please make sure the Mapping Toolbox is installed',...
                    ' so you can use polyxpoly.'])
            end
        end
        
        %% setup % reset function just reset time, setup function should reset everything
        %         function reset(W,seed)
        %             if ~exist('seed','var')
        %                 W.setup();
        %                 W.reset@world();
        %             else
        %                 rng(seed)
        %                 W.setup();
        %                 W.reset@world();
        %             end
        %         end
        function setup(W,seed)
            if exist('seed','var')
                % to see progress
                rng(seed)
            end
            W.placeCars();
            %             W.envCars = envCars;
            W.obstacle_history = [];
            % get room bounds
            B = W.bounds ;
            xlo = B(1) ; xhi = B(2) ; ylo = B(3) ; yhi = B(4) ;
            W.bounds_as_polyline = make_box([xhi-xlo,yhi-ylo]) + repmat(mean([xhi xlo ; yhi ylo],2),1,5) ;
            
            % get obstacle info
            %                 obs_size = W.obstacle_size_bounds ;
            %                 obs_rotation_bounds = W.obstacle_rotation_bounds ;
            
            % generate start position on left side of room with initial
            % heading of 0, and make sure it's not too close to the walls
            b = W.buffer ;
            
            xlo = xlo + 2*b ;
            xhi = xhi - 2*b ;
            ylo = ylo + 2*b ;
            yhi = yhi - 2*b ;
            
            if isempty(W.start)
                s = [xlo ;
                    rand_range(ylo, yhi) ;
                    0 ] ;
                W.start = s ;
            end
            
            % generate goal position on right side of room
            if isempty(W.goal)
                g = [xhi - 10 ;
                    (ylo+yhi)/2] ;
                W.goal = g ;
            end
            if ~isempty(W.envCars)
                W.N_obstacles = size(W.envCars,1)-1;
                % generate obstacles around room
                N_obs = W.N_obstacles ;
                
                if N_obs > 0
                    
                    O = nan(2, 6*N_obs) ; % preallocate obstacle matrix
                    
                    %                     W.obs_FRS =  cell(30,1);
                    %                     for i = 1:30
                    %                         W.obs_FRS {i} = nan(2, 6*N_obs) ;
                    %                     end
                    
                    %                     llo = obs_size(1) ; lhi = obs_size(2) ;
                    %                     orlo = obs_rotation_bounds(1) ;
                    %                     orhi = obs_rotation_bounds(2) ;
                    
                    %                     xlo = B(1) ; xhi = B(2) ; ylo = B(3) ; yhi = B(4) ;
                    %                     xlo = xlo + b ; xhi = xhi - b ;
                    %                     ylo = ylo + b ; yhi = yhi - b ;
                    
                    % obstacle rotation
                    % obstacle base
                    l = W.obstacle_size(1);
                    w = W.obstacle_size(2);
                    
                    o = [-l/2  l/2 l/2 -l/2 -l/2 ;
                        -w/2 -w/2 w/2  w/2 -w/2 ] ;
                    W.o_no_frs = o;
                    W.o = cell(length(W.v_array),30); % 10 is time index
                    %                     for v_idx = 1:length(W.v_array)
                    %                         for t_idx = 1: 30 % 1 second
                    %                             d = t_idx*0.1*W.v_array(v_idx);
                    %                             lo = 1/2*(2)*(t_idx*0.1)^2; %1/2*a*t^2;
                    %                             ro= 1/2*(2)*(t_idx*0.1)^2;
                    %                             W.o{v_idx,t_idx} = o+[d-lo d+ro d+ro d-lo d-lo; 0 0 0 0 0];
                    %                         end
                    %                     end
                    %
                    obs_count = 2; %start with second one, first one is ego
                    for idx = 1:6:(6*N_obs-1)
                        c = [W.envCars(obs_count,1); W.envCars(obs_count,3)];
                        O(:,idx:idx+4) = W.o_no_frs + repmat(c,1,5) ;
                        %                         [~,v_idx] = min(abs(W.v_array-W.envCars(obs_count,2)));
                        %                         for t_idx = 1:30
                        %                             W.obs_FRS{t_idx}(:,idx:idx+4) = W.o {v_idx,t_idx} + repmat(c,1,5) ;
                        %                         end
                        obs_count = obs_count + 1;
                    end
                    %
                    %
                    %                     if isnan(O(1,end))
                    %                         O = O(:,1:end-1) ;
                    %                     end
                    
                    W.obstacles = O ;
                    W.obstacles_old = O ;
                    
                    %                 W.obstacles_seen = [] ;
                    W.N_obstacles = N_obs ;
                end
                
                W.obstacles_unseen = [] ;
                W.obstacles_seen = W.obstacles;
                %
                %                 % set up plot data
                W.plot_data.obstacles_seen = W.obstacles ;
                W.plot_data.obstacles_unseen = [] ;
                W.plot_data.start = [] ;
                W.plot_data.goal = [] ;
                W.plot_data.goal_zone = [] ;
                W.plot_data.bounds = [] ;
                W.reset();
            end
        end
        
        function placefixedCars(W)
         W.envCars = zeros(W.num_cars, 6);
         xPos = [-85 -55 -20]
        laneNum = [2,1,2];
        for i = 1:W.num_cars
            W.envCars(i,1) = 100;
        W.envCars(i,2) = 0;%S.car_min_spd+rand*(S.car_max_spd-S.car_min_spd);
        W.envCars(i,3) =  2;
        W.envCars(i,4) = 0;%S.car_max_spd-rand*(S.car_max_spd-S.car_min_spd)/3;
        W.envCars(i,5) = 0
        W.envCars(i,6) = W.A_MM;%act #7 i
        end
        for i = 1:3
        W.envCars(end-i,1) = xPos(i);
        W.envCars(end-i,2) = 0;%S.car_min_spd+rand*(S.car_max_spd-S.car_min_spd);
        W.envCars(end-i,3) = 4*laneNum(i) + 2;
        W.envCars(end-i,4) = 0;%S.car_max_spd-rand*(S.car_max_spd-S.car_min_spd)/3;
        W.envCars(end-i,5) = laneNum(i);
        W.envCars(end-i,6) = W.A_MM;%act #7 is maintain, for env cars only some actions are available
        end     
        
        end
        function placeCars(W)
            W.envCars = zeros(W.num_cars, 6);
            
            MAX_DIST =  W.SIM_MAX_DISTANCE - W.car_safe_dist;
            for i = 2:W.num_cars
                laneOverlap = true;
                while laneOverlap
                    xPos = -MAX_DIST * rand  + W.car_safe_dist;
                    laneNum = randi([0,2]);
                    laneOverlap = false;
                    for j = 1:i
                        %if overlap, y dir close and x dir close
                        % TODO: check with eco car as well, ego not here
                        if  abs(xPos - W.envCars(j,1)) < W.car_safe_dist %not useful any more, static obstacles%&&laneNum ==S.envCars(j,5)
                            laneOverlap = true;
                            break;
                        end
                    end
                    if ~laneOverlap % assign car to new slot
                        W.envCars(i,1) = xPos;
                        W.envCars(i,2) = 0;%S.car_min_spd+rand*(S.car_max_spd-S.car_min_spd);
                        W.envCars(i,3) = 4*laneNum + 2;
                        W.envCars(i,4) = 0;%S.car_max_spd-rand*(S.car_max_spd-S.car_min_spd)/3;
                        W.envCars(i,5) = laneNum;
                        W.envCars(i,6) = W.A_MM;%act #7 is maintain, for env cars only some actions are available
                        %action = 1 2 3 4 5  6 7 8 9  10 11 12;
                        %         L L L L M  M M M R   R  R  R
                        %        2d d 0 a 2d d m a 2d  d  m  a
                    end
                end
            end
        end
        
        function update_envCars(W,envCars)
            
            W.envCars = envCars;
            
            
            W.obstacles_old = W.obstacles;
            %             obs_count = 2;
            %             for idx = 1:6:(6*W.N_obstacles-1)
            %                 c = [W.envCars(obs_count,1); W.envCars(obs_count,3)];
            %                 W.obstacles(:,idx:idx+4) = W.o_no_frs + repmat(c,1,5) ;
            %                 [~,v_idx] = min(abs(W.v_array-W.envCars(obs_count,2)));
            %                 for t_idx = 1:30
            %                     W.obs_FRS{t_idx}(:,idx:idx+4) = W.o {v_idx,t_idx} + repmat(c,1,5) ;
            %                 end
            %                 obs_count = obs_count + 1;
            %             end
            %             W.obstacles_seen = W.obstacles;
        end
        
        function world_info = get_world_info(W,agent_info,~)
            W.vdisp('Getting world info!',3)
            
            if nargin > 1
                zcur = agent_info.state(agent_info.position_indices,end)' ;
                zcur = round(zcur,6) ;
                z = unique(zcur,'rows')' ;
                r = agent_info.sensor_radius ;
            else
                W.vdisp('No agent info provided!',2)
                z = W.start ;
                r = 1 ;
            end
            
            %             O = W.obstacles_unseen ;
            %
            %             if ~isempty(O)
            %                 N = size(O,2)/6 ; % number of obstacles still unseen
            %
            %                 indices_seen = [] ;
            %
            %                 O_out = [] ;
            %                 for idx = 1:6:(6*N-1)
            %                     if size(z,2) == 1
            %                         dToObs = dist_point_to_polyline(z,O(:,idx:idx+4)) ;
            %                     else
            %                         dToObs = poly_poly_dist(z(1,:),z(2,:),...
            %                             O(1,idx:idx+4),O(2,idx:idx+4)) ;
            %                     end
            %                     if any(dToObs <= r)
            %                         O_out = [O_out, O(:,idx:idx+4),nan(2,1)] ;
            %                         indices_seen = [indices_seen, idx:idx+5] ;
            %                     end
            %
            %                     % make sure the last index seen is within the length of O
            %                     if ~isempty(indices_seen) && indices_seen(end) > size(O,2)
            %                         indices_seen = indices_seen(1:end-1) ;
            %                     end
            %                 end
            
            % if some obstacles were seen, remove them from O
            %                 O(:,indices_seen) = [] ;
            %                 W.obstacles_unseen = O ;
            %                 W.obstacles_seen = [W.obstacles_seen, O_out] ;
            %             end
            
            world_info.obstacles = W.obstacles  ;
            world_info.bounds = W.bounds ;
            world_info.start = W.start ;
            world_info.goal = W.goal ;
            world_info.dimension = W.dimension ;
            world_info.obs_FRS = W.obs_FRS;
        end
        
        
        %% collision check
        
        function [signed_dist, rel_spd, rel_y, ego_vy] = getCar_abs(W,i,agent_info)
            
            cars = W.envCars(2:end,:);
            ego_x  = agent_info.state(1,end);
            ego_v  = agent_info.state(4,end);
            ego_y  = agent_info.state(2,end);
            try
                ego_vy = (agent_info.state(2,end)-agent_info.state(2,end-1))/(agent_info.time(end) - agent_info.time(end-1)); %ynew-yold / dT
            catch E
                warning('Undefined ego_vy!')
                ego_vy = 0;
            end
            if i == 1
                lane_cars = cars(cars(:,5) == 2,:);
                cars = lane_cars(lane_cars(:,1) < ego_x,:);
                cary = 10;
            elseif i == 2
                lane_cars = cars(cars(:,5) ==  2,:);
                cars = lane_cars(lane_cars(:,1) > ego_x,:);
                cary = 10;
            elseif i == 3
                lane_cars = cars(cars(:,5) == 1,:);
                cars = lane_cars(lane_cars(:,1) < ego_x,:);
                cary = 6;
            elseif i == 4
                lane_cars = cars(cars(:,5) == 1,:);
                cars = lane_cars(lane_cars(:,1) > ego_x,:);
                cary = 6;
            elseif i == 5
                lane_cars = cars(cars(:,5) == 0,:);
                cars = lane_cars(lane_cars(:,1) < ego_x,:);
                cary = 2;
            elseif i == 6
                lane_cars = cars(cars(:,5) == 0,:);
                cars = lane_cars(lane_cars(:,1) > ego_x,:);
                cary = 2;
            end
            [~,idx] = min(abs(cars(:,1) - ego_x));
            wanted_car = cars(idx,:);
            %1   2
            %3 e 4
            %5   6
            if ~isempty(wanted_car) && abs(wanted_car(1)-ego_x) < W.car_max_vis_dist
                signed_dist = wanted_car(1)-ego_x;
                rel_spd = wanted_car(2)-ego_v;
            else
                if i == 2 || i ==4 || i == 6
                    signed_dist = W.car_max_vis_dist;
                    rel_spd = W.car_max_spd - ego_v;
                else
                    signed_dist = -W.car_max_vis_dist;
                    rel_spd = W.car_min_spd - ego_v;
                end
            end
            rel_y     =  cary - ego_y;
        end
        function [total_rew]= getRew(W,agent_info,observation)
            %observation(2)
            ego_y = agent_info.state(2, end);
            diff_y = observation(2);
            y_near_obst = diff_y + ego_y;
            diff_y = abs(diff_y);
            if (diff_y <= 10) && (diff_y >= 2)
                Reward_Y = 5 * exp(- 1/(diff_y.^2+1));
            elseif (diff_y < 2)
                Reward_Y = 5 * exp(- 1/((y_near_obst-2)^2+1)) - 2 * abs(abs(diff_y)-2) ;
            elseif (diff_y > 10)
                Reward_Y = 5 * exp(- 1/((y_near_obst-10)^2+1)) - 2 * abs(abs(diff_y)-10) ;
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
            Reward_speed = 4 * exp(- 1./(observation(5)^2+1)) - 3.7;
%             total_rew = Reward_speed;
            total_rew = (0.8*Reward_speed+Reward_Y-4);  %reward scaled by speed
        end
        function Observation = get_ob_deprecated(W,agent_info)
            Observation = zeros(5,1);
            min_index = 2;
            min_dis = inf;
            min_y = 0;
            
            for idx_neighbor = [2,4,6] % observation here can do better, use a hybrid of desired pose and actual pose;
                %currently using desired pose since it will eventually converge to desired pose and current pose doesnt really matter
                [rel_dist,rel_v,rel_y, ego_vy]=W.getCar_abs(idx_neighbor,agent_info);
                if (rel_dist < min_dis)
                    min_index = idx_neighbor;
                    min_dis = rel_dist;
                    min_y = rel_y;
                end
            end
            Observation(1) = min_dis;
            Observation(2) = min_y;
            
            %             iniOb(19:21) = [S.envCars(1,2);S.envCars(1,3);ego_vy];
            %iniOb(1) = S.W.goal(1) -  S.envCars(1,1);
            Observation(3) = agent_info.state(1,end); %ego_x
            Observation(4) = agent_info.state(2,end); %ego_y
            Observation(5) = agent_info.state(4,end); %ego_vx
            
        end
        function Observation = get_ob(W,agent_info)
            Observation = zeros(6,1);
            rel_dist_vec = [];
            rel_y_vec = [];
            for idx_neighbor = [2,4,6] % observation here can do better, use a hybrid of desired pose and actual pose;
                %currently using desired pose since it will eventually converge to desired pose and current pose doesnt really matter
                [rel_dist,~,rel_y,~]=W.getCar_abs(idx_neighbor,agent_info);
                rel_dist_vec = [rel_dist_vec;rel_dist];
                rel_y_vec = [rel_y_vec;rel_y];
%                 if (rel_dist < min_dis)
% %                     min_index = idx_neighbor;
%                     min_dis = rel_dist;
%                     min_y = rel_y;
%                 end
            end
            [min_dis, min_idx] = min(rel_dist_vec);
            min_rel_y = rel_y_vec(min_idx);
            rel_dist_vec(min_idx) = [];
            rel_y_vec(min_idx) = [];
            [second_min_dis, second_min_idx] = min(rel_dist_vec);
            second_min_rel_y = rel_y_vec(second_min_idx);

            Observation(1) = min_dis;
            Observation(2) = min_rel_y;
            Observation(3) = second_min_dis - min_dis;
            Observation(4) = second_min_rel_y;
            %             iniOb(19:21) = [S.envCars(1,2);S.envCars(1,3);ego_vy];
            %iniOb(1) = S.W.goal(1) -  S.envCars(1,1);
%             Observation(3) = agent_info.state(1,end); %ego_x
%             Observation(4) = agent_info.state(2,end); %ego_y
            Observation(5) = agent_info.state(4,end); %ego_vx
            Observation(6) = agent_info.state(2,end); %y
            
        end
        
        function out = collision_check(W,agent_info,check_full_traj_flag)
            
            % by default, don't check the full trajectory
            if nargin < 3
                check_full_traj_flag = false ;
            end
            
            % initialize output (innocent until proven guilty)
            out = 0 ;
            
            % extract agent info
            pos_idx = agent_info.position_indices ;
            h_idx = agent_info.heading_index ;
            fp = agent_info.footprint_vertices ;
            Z = agent_info.state ;
            T = agent_info.time ;
            
            % set up obstacles
            O = [W.obstacles, nan(2,1), W.bounds_as_polyline,nan(2,1)] ;
            O_old = [W.obstacles_old, nan(2,1), W.bounds_as_polyline,nan(2,1)] ;
            % if the agent has moved, we need to check if it crashed;
            % alternatively, we could check the entire trajectory for crashes;
            % finally, if it didn't move, then check the full trajectory
            t_start = W.current_time ;
            if check_full_traj_flag
                t_log = true(size(T)) ;
            else
                t_log = T >= t_start ;
            end
            
            if sum(t_log) > 1
                T = T(t_log) ;
                Z = Z(:,t_log) ;
                
                % create the desired time vector; check for collision every
                % 10ms by default, in 5s chunks of time
                dt = T(end) - T(1) ; % total time of traj
                N_chk = ceil(dt/5) ; % number of checks to make
                
                t_world = T(1) ; % start time of each check
                
                for chk_idx = 1:N_chk
                    % get the time vector to check
                    t_log_idx = T >= t_world ;
                    T_idx = T(t_log_idx) ;
                    Z_idx = Z(:,t_log_idx) ;
                    
                    % create the time vector for interpolation
                    t_chk = (0:0.02:5 ) + t_world ;
                    t_chk = t_chk(t_chk <= T_idx(end)) ;
                    
                    % create the interpolated trajectory
                    %                     try
                    Z_chk = match_trajectories(t_chk,T_idx,Z_idx) ;
                    %                         if ~empty(W.obstacles_old)
                    O_idx = 1;
                    while O_idx + 5 <= size(O_old,2)
                        if abs(O(1,O_idx)-O_old(1,O_idx)) > 100
                            O(:,O_idx:O_idx+5) = []; O_old(:,O_idx:O_idx+5) = [];
                        else
                            O_idx = O_idx + 6;
                        end
                    end
                    O_x = match_trajectories(t_chk,[T_idx(1) T_idx(end)],[O_old(1,:)',O(1,:)'])' ;
                    O_y = match_trajectories(t_chk,[T_idx(1) T_idx(end)],[O_old(2,:)',O(2,:)'])' ;
                    
                    %                         end
                    
                    % create a copy of the agent footprint, rotated at each
                    % point of the trajectory
                    X = Z_chk(pos_idx,:) ; % xy positions of trajectory
                    N = size(X,2) ;
                    X = repmat(X(:),1,size(fp,2)) ;
                    
                    % X is rep'd for each pt of agent footprint vertices
                    
                    if ~isempty(h_idx) && (length(agent_info.footprint) > 1)
                        % if there is a heading, and the agent is not a circle,
                        % then rotate each footprint contour
                        H = Z_chk(h_idx,:) ;
                        R = rotmat(double(H)) ;
                        F = R*repmat(fp,N,1) + X ;
                    else
                        % otherwise, just place the footprint contour at each
                        % point of the trajectory
                        F = repmat(fp,N,1) + X ;
                    end
                    
                    Fx = [F(1:2:end,:)' ; nan(1,N)] ;
                    Fy = [F(2:2:end,:)' ; nan(1,N)] ;
                    F = [Fx(:)' ; Fy(:)'] ;
                    ci_bool = ones( 1,size(F,2)/6);
                    for chk_time_idx = 1: size(F,2)/6
                        F_idx = (6*(chk_time_idx-1)+1):(6*(chk_time_idx));
                        [ci,cyi] =polyxpoly(F(1,F_idx),F(2,F_idx),O_x(chk_time_idx,:),O_y(chk_time_idx,:)) ;
                        %                         figure(4); clf
                        %                         mapshow(F(1,F_idx),F(2,F_idx),'DisplayType','polygon','LineStyle',':')
                        %                         mapshow(O_x(chk_time_idx,:),O_y(chk_time_idx,:),'Marker','+')
                        %                         xlim([F(1,1)-50 F(1,1)+100]);                          mapshow(ci,cyi,'DisplayTcype','point','Marker','o')
                        if ~isempty(ci)
                            
                            break
                        else
                            ci_bool(chk_time_idx)   = 0;
                        end
                    end
                    
                    if any(ci_bool)
                        out = 1 ;
                    else
                        out = 0;
                    end
                    % check if the resulting contour intersects the obstacles
                    %                         [ci,cyi] = polyxpoly(F(1,:),F(2,:),O(1,:),O(2,:)) ;
                    
                    
                    
                    
                    % increment the time index
                    t_world = t_world + 5 ;
                    %                     catch
                    %                         W.vdisp('Check failed, skipping to next portion!',2)
                    %                         out = -1 ;
                    %                     end
                end
            end
            
            % update the world time index
            W.current_time = agent_info.time(end) ;
        end
        
        %% plotting
        function plot(W,~)
            % set up hold if needed
            hold_check = false ;
            if ~ishold
                hold_check = true ;
                hold on
            end
            
            % plot sensed obstacles
            O_seen = W.obstacles_seen ;
            
            if isempty(O_seen)
                O_seen = nan(2,1) ;
            end
            
            if ~check_if_plot_is_available(W,'road_lanes')
                road_lanes ={};
                w= 2;
                road_lanes{1} =  fill([-W.start_line 0 0 -W.start_line -W.start_line],[-0.7 -0.7 12.7 12.7 -0.7],[207,207,207]/255);
                % plot([0,15],[2, 2],'--w','LineWidth',8)
                road_lanes{2} = plot([-W.start_line,0],[12, 12],'LineWidth',w,'Color',[255, 255, 255]/255);
                road_lanes{3} =  plot([-W.start_line,0],[8, 8],'--','LineWidth',w,'Color',[1 1 1]);
                road_lanes{4} =  plot([-W.start_line,0],[4, 4],'--','LineWidth',w,'Color',[1 1 1]);
                road_lanes{5} =  plot([-W.start_line,0],[0, 0],'LineWidth',w,'Color',[255, 255, 255]/255);
                W.plot_data.road_lanes = road_lanes ;
            end
            
%             if check_if_plot_is_available(W,'obstacles_seen')
%                 W.plot_data.obstacles_seen.XData = O_seen(1,:) ;
%                 W.plot_data.obstacles_seen.YData = O_seen(2,:) ;
%             else
           % replot obstalces every time.
           [X_plot, Y_plot] = W.convert_data_for_plots([nan O_seen(1,:)],[nan O_seen(2,:)]);
           
           ph=patch(X_plot, Y_plot, W.obstacle_seen_color);
           ph.FaceAlpha = W.obstacle_alpha;
%                 seen_data = plot(O_seen(1,:),O_seen(2,:),'Color',W.obstacle_seen_color) ;
%                 W.plot_data.obstacles_seen = seen_data ;
%             end
            
            % plot unsensed obstacles
            O_unseen = W.obstacles_unseen ;
            
            if isempty(O_unseen)
                O_unseen = nan(2,1) ;
            end
            
            if check_if_plot_is_available(W,'obstacles_unseen')
                W.plot_data.obstacles_unseen.XData = O_unseen(1,:) ;
                W.plot_data.obstacles_unseen.YData = O_unseen(2,:) ;
            else
                unseen_data = plot(O_unseen(1,:),O_unseen(2,:),'Color',W.obstacle_unseen_color) ;
                W.plot_data.obstacles_unseen = unseen_data ;
            end
            
            % plot start
            s = W.start ;
            
%             if ~check_if_plot_is_available(W,'start')
%                 start_data = plot(s(1),s(2),'bx') ;
%                 W.plot_data.start = start_data ;
%             end
            
            % plot goal and goal zone
            g = W.goal ;
            g_circ = make_circle(W.goal_radius,100) + repmat(g(1:2),1,100) ;
            
            if ~check_if_plot_is_available(W,'goal')
                goal_data = plot(g(1),g(2),'x','Color',[0.2 0.7 0]) ;
                W.plot_data.goal = goal_data ;
            end
            
            if ~check_if_plot_is_available(W,'goal_zone')
                goal_zone_data = plot(g_circ(1,:),g_circ(2,:),'--','Color',[0.2 0.7 0]) ;
                W.plot_data.goal_zone = goal_zone_data ;
            end
            
            % plot bounds
            B = W.bounds_as_polyline ;
            
            if ~check_if_plot_is_available(W,'bounds')
                bounds_data = plot(B(1,:),B(2,:),'Color',W.obstacle_seen_color) ;
                W.plot_data.bounds = bounds_data ;
            end
            
            if hold_check
                hold off ;
            end
        end
        
        function plot_at_time(W,~)
            % set up hold if needed
            hold_check = false ;
            if ~ishold
                hold_check = true ;
                hold on
            end
            
            % plot ass obstacles
            O = W.obstacles ;
            
            if isempty(O)
                O = nan(2,1) ;
            end
            
            if check_if_plot_is_available(W,'obstacles')
                W.plot_data.obstacles.XData = O(1,:) ;
                W.plot_data.obstacles.YData = O(2,:) ;
            else
                seen_data = plot(O(1,:),O(2,:),'Color',W.obstacle_seen_color) ;
                W.plot_data.obstacles = seen_data ;
            end
            
            % plot start
            s = W.start ;
            
            if ~check_if_plot_is_available(W,'start')
                start_data = plot(s(1),s(2),'bx') ;
                W.plot_data.start = start_data ;
            end
            
            % plot goal and goal zone
            g = W.goal ;
            g_circ = make_circle(W.goal_radius,100) + repmat(g(1:2),1,100) ;
            
            if ~check_if_plot_is_available(W,'goal')
                goal_data = plot(g(1),g(2),'x','Color',[0.2 0.7 0]) ;
                W.plot_data.goal = goal_data ;
            end
            
            if ~check_if_plot_is_available(W,'goal_zone')
                goal_zone_data = plot(g_circ(1,:),g_circ(2,:),'--','Color',[0.2 0.7 0]) ;
                W.plot_data.goal_zone = goal_zone_data ;
            end
            
            % plot bounds
            B = W.bounds_as_polyline ;
            
            if ~check_if_plot_is_available(W,'bounds')
                bounds_data = plot(B(1,:),B(2,:),'Color',W.obstacle_seen_color) ;
                W.plot_data.bounds = bounds_data ;
            end
            
            if hold_check
                hold off ;
            end
            
            % set axes
            axis equal
        end
        function [X, Y] = convert_data_for_plots(W,x,y)
            idn=find(isnan(x));
            Sz=diff(idn)-1;
            Nmax=max(Sz);
            N=numel(Sz);
            X=zeros(Nmax,N);
            Y=X;
            for i=1:N
                X(1:Sz(i),i)=x(idn(i)+1:idn(i+1)-1);
                Y(1:Sz(i),i)=y(idn(i)+1:idn(i+1)-1);
            end
        end
    end
    
end