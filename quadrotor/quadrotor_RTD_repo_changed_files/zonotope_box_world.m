classdef zonotope_box_world < world
properties
    % goal info
    goal_zone = [] ; % this is a struct used to plot the goal
    goal_type = 'sphere' ; % 'sphere' or 'box'
    
    % start info
    start_zone = [] ; % this is a struct used to plot the start

    % obstacle info
    use_wall_obstacles_flag = false ;
    buffer_start = 5 ;
    buffer_goal = 5 ;
    N_tall = 15 ;
    N_wide = 10 ;
    N_long = 5 ;
    N_boxy = 15 ;
    obs_tall_size_range = [0.25, 0.5] ; % [min,max] thickness
    obs_wide_size_range = [0.25, 0.5] ; % [min,max] thickness
    obs_long_size_range = [0.25, 0.5] ; % [min,max] thickness
    obs_boxy_size_range = [2.0, 3.0, 2.5, 1] ; % [min,max,mean,std] side length
    create_obstacle_center_timeout = 5 ; % s
    

    % plot
    plot_obs_face_color = [1 0.7 0.7] ;
    plot_obs_face_opacity = 0.7 ;
    plot_obs_sensed_face_color = [1 0.6 0.6] ;
    
    plot_obs_edge_color = [0.5 0.1 0.1] ;
    plot_obs_edge_opacity = 0.7 ;
    plot_obs_sensed_edge_color = [0.5 0 0] ;
    plot_sensed_obstacles_flag = true ;
    
    plot_ground_color = [0.7 0.5 0.4] ;
    plot_ground_opacity = 0.1 ;
    plot_ceiling_color = [1 0.8 0.8] ;
    plot_ceiling_opacity =0;
    plot_wall_color = [1 0.8 0.8] ;
    plot_wall_opacity = 0;
    
    plot_goal_zone_color = [0 1 0.5] ;
    plot_goal_zone_opacity = 0.1 ;
    plot_start_zone_color = [0.6 0 0.7] ;
    plot_start_zone_opacity = 0.0 ;
    
    plot_goal_zone_flag = true ;
    plot_start_zone_flag = true ;
end
methods
    %% constructor
    function W = zonotope_box_world(varargin)
        % default properties
%         W.bounds = [0 25 -2.5 2.5 -2.5 2.5] ;
        W.bounds = [0 100 -2.5 2.5 -2.5 2.5] ;
        W.dimension = 3 ;
        W.start = [1;
            0.5*(W.bounds(4) + W.bounds(3)) ;
            0.5*(W.bounds(6) + W.bounds(5))] ;
        W.goal_radius = 2.5 ;
        W.obstacles = {} ;

        % parse arguments
        W = parse_args(W,varargin{:}) ;
        W.start = W.start(:) ;
        W.goal = W.goal(:) ;
        
        % make goal
        if length(W.goal) ~= 3
            r = W.goal_radius ;
            W.goal = [W.bounds(2) - r;
                rand_range(W.bounds(3)+r,W.bounds(4)-r);
                rand_range(W.bounds(5)+r,W.bounds(6)-r)] ;
        end
        
        % set up for plotting
        W.make_start_zone() ;
        W.make_goal_zone() ;
        W.plot_data.start = [] ;
        W.plot_data.start_zone = [] ;
        W.plot_data.goal = [] ;
        W.plot_data.goal_radius = [] ;
        W.plot_data.goal_zone = [] ;
        
        % make obstacles
        if isempty(W.obstacles)
            W.setup() ;
        else
            W.N_obstacles = length(W.obstacles) ;
        end
    end

    %% setup
    function setup(W,seed)
        if exist('seed','var')
            rng(seed)
        end
        %reset obstacles
        W.obstacles = {};
        
        % create ground and ceiling obstacles
        W.obstacles{1} = W.create_ground() ;
        W.obstacles{2} = W.create_ceiling() ;
        
        % create wall obstacles
        if W.use_wall_obstacles_flag
            W.create_walls() ;
        end

        idx = length(W.obstacles) + 1 ; % idx to keep track of total
        
        % create tall obstacles
        for oidx = 1:W.N_tall
            obs = [] ;
            while isempty(obs)
                obs = W.create_obstacle('tall') ;
            end
            W.obstacles{idx} = obs ;
            idx = idx + 1 ;
        end

        % create wide obstacles
        for oidx = 1:W.N_wide
            obs = [] ;
            while isempty(obs)
                obs = W.create_obstacle('wide') ;
            end
            W.obstacles{idx} = obs ;
            idx = idx + 1 ;
        end

        % create long obstacles
        for oidx = 1:W.N_long
            obs = [] ;
            while isempty(obs)
                obs = W.create_obstacle('long') ;
            end
            W.obstacles{idx} = obs ;
            idx = idx + 1 ;
        end

        % create boxy obstacles
        for oidx = 1:W.N_boxy
            obs = [] ;
            while isempty(obs)
                obs = W.create_obstacle('boxy') ;
            end
            W.obstacles{idx} = obs ;
            idx = idx + 1 ;
        end

        W.N_obstacles = length(W.obstacles) ;

        W.reset()
    end
    
    function reset(W)
        % reset time
        W.current_time = 0 ;

        % reset plot data
        W.plot_data.start = [] ;
        W.plot_data.start_zone = [] ;
        W.plot_data.goal = [] ;
        W.plot_data.goal_radius = [] ;
        W.plot_data.goal_zone = [] ;
    end

    %% create obstacles
    function add_obstacle(W,varargin)
        % W.add_obstacle(obs_type)
        % W.add_obstacle(length,width,height,center)
        %
        % Adds a static obstacle either of the type passed in (tall, wide,
        % long, or boxy),
        
        N_in = length(varargin) ;
        if N_in == 1
            obs = W.create_obstacle(varargin{1}) ;
        else
            l = varargin{1} ;
            w = varargin{2} ;
            h = varargin{3} ;
            if N_in == 4
                c = varargin{4} ;
            else
                c = zeros(3,1) ;
            end
            [l,w,h,c] = resize_box_for_world_bounds(l,w,h,c,W.bounds) ;
            obs = zonotope_obstacle(l,w,h,c) ;
        end
        
        W.update_obstacle_list(obs) ;
    end
    
    function update_obstacle_list(W,obs)
        W.obstacles = [W.obstacles, {obs}] ;
        W.N_obstacles = W.N_obstacles + 1 ;
    end
    
    function obs = create_obstacle(W,obs_type)
        switch obs_type
            case 'tall'
                [l,w,h,c] = W.create_tall_obstacle() ;
            case 'wide'
                [l,w,h,c] = W.create_wide_obstacle() ;
            case 'long'
                [l,w,h,c] = W.create_long_obstacle() ;
            case 'boxy'
                [l,w,h,c] = W.create_boxy_obstacle() ;
            otherwise
                error('Obstacle type must be tall, wide, long, or boxy!')
        end
        
        B_adjusted = W.bounds ;
        B_adjusted(1) = B_adjusted(1) + W.buffer_start ;
        B_adjusted(2) = B_adjusted(2) - W.buffer_goal ;
        [l,w,h,c] = resize_box_for_world_bounds(l,w,h,c,B_adjusted) ;

        if ~isempty(l)
            obs = zonotope_obstacle(l,w,h,c) ;

            obs.plot_face_color = W.plot_obs_face_color ;
            obs.plot_face_opacity = W.plot_obs_face_opacity ;

            obs.plot_edge_color = W.plot_obs_edge_color ;
            obs.plot_edge_opacity = W.plot_obs_edge_opacity ;
        else
            warning('Obstacle creation failed!')
            obs = [] ;
        end
    end

    function [l,w,h,c] = create_tall_obstacle(W)
        W.vdisp('Creating tall obstacle',7)

        dlo = W.obs_tall_size_range(1) ;
        dhi = W.obs_tall_size_range(2) ;

        l = rand_range(dlo,dhi) ;
        w = rand_range(dlo,dhi) ;
        h = W.bounds(6) - W.bounds(5) ;
        
        c = W.create_obstacle_center(l,w,h,true) ;
    end

    function [l,w,h,c] = create_wide_obstacle(W)
        W.vdisp('Creating wide obstacle',7)

        dlo = W.obs_wide_size_range(1) ;
        dhi = W.obs_wide_size_range(2) ;

        l = rand_range(dlo,dhi) ;
        w = W.bounds(4) - W.bounds(3) ;
        h = rand_range(dlo,dhi) ;
        
        c = W.create_obstacle_center(l,w,h,false) ;
    end

    function [l,w,h,c] = create_long_obstacle(W)
        W.vdisp('Creating long obstacle',7)

        dlo = W.obs_long_size_range(1) ;
        dhi = W.obs_long_size_range(2) ;

        l = W.bounds(2) - W.bounds(1) ;
        w = rand_range(dlo,dhi) ;
        h = rand_range(dlo,dhi) ;
        
        c = W.create_obstacle_center(l,w,h,false) ;
    end

    function [l,w,h,c] = create_boxy_obstacle(W)
        W.vdisp('Creating boxy obstacle',7)

        dlo = W.obs_boxy_size_range(1) ;
        dhi = W.obs_boxy_size_range(2) ;
        m = W.obs_boxy_size_range(3) ;
        s = W.obs_boxy_size_range(4) ;

        l = rand_range(dlo,dhi,m,s) ;
        w = rand_range(dlo,dhi,m,s) ;
        h = rand_range(dlo,dhi,m,s) ;
        
        c = W.create_obstacle_center(l,w,h,false) ;
    end

    function center = create_obstacle_center(W,l,w,h,center_vertically_flag)
        if nargin < 5
            center_vertically_flag = true ;
        end

        t_start = tic ;
        center = zeros(3,1) ;
        
        center(1) = rand_range((W.bounds(1) + l/2 + W.buffer_start),...
                               (W.bounds(2) - l/2 - W.buffer_goal)) ;
        center(2) = rand_range(W.bounds(3), W.bounds(4)) ;
        center(3) = rand_range(W.bounds(5), W.bounds(6)) ;
        
        if center_vertically_flag
            center(3) = (W.bounds(6) + W.bounds(5))/2 ;
        end
        
        d_to_start = dist_point_to_box(W.start,l,w,h,center) ;
        d_to_goal = dist_point_to_box(W.goal,l,w,h,center) ;
        
        while (d_to_start < W.buffer_start) || (d_to_goal < (W.buffer_goal + W.goal_radius))
            center(1) = rand_range(W.bounds(1), W.bounds(2)) ;
            center(2) = rand_range(W.bounds(3), W.bounds(4)) ;
            if ~center_vertically_flag
                center(3) = rand_range(W.bounds(5), W.bounds(6)) ;
            end
            
            d_to_start = dist_point_to_box(W.start,l,w,h,center) ;
            d_to_goal = dist_point_to_box(W.goal,l,w,h,center) ;
            
            if toc(t_start) > W.create_obstacle_center_timeout
                error('Timed out while trying to create obstacle center!')
            end
        end
    end
    
    function obs = create_ground(W)
        % create ground obstacles
        obs = zonotope_obstacle(W.bounds(2) - W.bounds(1),...
                                            W.bounds(4) - W.bounds(3),...
                                            0.1,...
                                            [(W.bounds(2) + W.bounds(1))/2;
                                            (W.bounds(4) + W.bounds(3))/2;
                                            W.bounds(5) + 0.05]);
        obs.plot_face_opacity = W.plot_ground_opacity ;
        obs.plot_face_color = W.plot_ground_color ;
    end
    
    function obs = create_ceiling(W)
        % create ceiling obstacle
        obs  = zonotope_obstacle(W.bounds(2) - W.bounds(1),...
                                            W.bounds(4) - W.bounds(3),...
                                            0.1,...
                                            [(W.bounds(2) + W.bounds(1))/2;
                                            (W.bounds(4) + W.bounds(3))/2;
                                            W.bounds(6) - 0.05]);
        obs.plot_face_opacity = W.plot_ceiling_opacity ;
        obs.plot_edge_opacity = W.plot_ceiling_opacity ;
        obs.plot_edge_color = W.plot_ceiling_color ;
    end
    
    function obs = create_wall(W,side)
        % obs = create_wall(W,side)
        %
        % Create a wall obstacle. The side input should be a string, either
        % 'left', 'right', 'front', or 'back'.
        
        B = W.bounds ;
        
        switch side
            case 'left'
                l = B(2) - B(1) ;
                w = 0.02 ;
                h = B(6) - B(5) ;
                c = [0.5*(B(2) + B(1)) ;
                     B(3) ;
                     0.5*(B(6) + B(5))] ;
                
            case 'right'
                l = B(2) - B(1) ;
                w = 0.02 ;
                h = B(6) - B(5) ;
                c = [0.5*(B(2) + B(1)) ;
                     B(4) ;
                     0.5*(B(6) + B(5))] ;
                
            case 'front'
                l = 0.02 ;
                w = B(4) - B(3) ;
                h = B(6) - B(5) ;
                c = [B(1) ;
                     0.5*(B(3) + B(4)) ;
                     0.5*(B(6) + B(5))] ;
                
            case 'back'
                l = 0.02 ;
                w = B(4) - B(3) ;
                h = B(6) - B(5) ;
                c = [B(2) ;
                     0.5*(B(3) + B(4)) ;
                     0.5*(B(6) + B(5))] ;
        end
        obs = zonotope_obstacle(l,w,h,c);
        
        obs.plot_face_opacity = W.plot_wall_opacity ;
        obs.plot_face_color = W.plot_wall_color ;
    end
    
    function create_walls(W)
        % create walls
        W.obstacles{3} = W.create_wall('left') ;
        W.obstacles{4} = W.create_wall('right') ;
        W.obstacles{5} = W.create_wall('front') ;
        W.obstacles{6} = W.create_wall('back') ;
    end

    %% get world info
    function world_info = get_world_info(W,agent_info,~)
        % get world info
        world_info.goal = W.goal ;
        world_info.goal_radius = W.goal_radius ;
        world_info.start = W.start ;
        world_info.bounds = W.bounds ;
        world_info.dimension = W.dimension ;
        
        % get agent current position
        x = agent_info.position(:,end) ;
        
        % return all obstacles that are within the sensor radius at the
        % current time
        O = {} ;
        oidx = 1 ;
        for idx = 1:W.N_obstacles
            o = W.obstacles{idx} ;
            
            % get the obstacle's vertices
            V = o.vertices ;
            
            % get the closest point on the obstacle to the current agent
            % position
            x_close = closest_point_on_box(x,V) ;
            
            % if the closest point is within the sensor radius, add that
            % obstacle to the list
            if vecnorm(x_close(:) - x(:)) <= agent_info.sensor_radius
                O{oidx} = o ;
                oidx = oidx + 1 ;
                
                if W.plot_sensed_obstacles_flag && idx > 2 % NOTE THIS IS A HACK FOR GROUND AND CEILING
                    W.obstacles{idx}.plot_face_color = W.plot_obs_sensed_face_color ;
                end
            end
        end
        
        % output
        world_info.obstacles = O ;
    end

    %% collision check
    function out = collision_check(W,agent_info,check_full_traj_flag)
        if nargin < 3
            check_full_traj_flag = false ;
        end

        % get agent trajectory
        [T,Z,R] = W.collision_check_get_agent_states(agent_info) ;
        V = agent_info.collision_check_input_data.vertices ;
        
        % get points to check for collision
        V_check = W.collision_check_point_prep(T,Z,R,V,check_full_traj_flag) ;
        
        if ~isempty(V_check)
            out = W.collision_check_loop(V_check) ;
        else
            W.vdisp('Agent time has not updated since last collision check.',3)
            out = false ;
        end

        % update world time
        W.current_time = agent_info.time(end) ;
    end
    
    function [T,Z,R] = collision_check_get_agent_states(W,agent_info)
        T = agent_info.time ;
        Z = agent_info.position ;
        R = agent_info.attitude ;
    end
    
    function V_check = collision_check_point_prep(W,T,Z,R,V,check_full_traj_flag)
        N_V = size(V,1) ; % number of vertices

        % get vertices into 3-by-N format, since they are usually going to be
        % passed in as N-by-3 (for plotting a patch)
        [r,c] = size(V) ;
        if c == 3 && r ~= 3 % otherwise it could be a triangle I guess
            V = V' ;
        end

        % get all the time indices that are relevant (approximately)
        if ~check_full_traj_flag
            t = W.current_time ;
            t_log = T >= t ;
        else
            t_log = true(size(T)) ;
        end

        N_check = sum(t_log) ;

        if N_check > 0
            Z = Z(:,t_log) ;
            R = R(:,:,t_log) ;
            T = T(t_log) ;

            % put vertices at each point along trajectory; first, get the rotation
            % matrices into a big matrix
            R_cel = mat2cell(R,3,3,ones(1,N_check)) ;
            R_blk = blkdiag(R_cel{:}) ;

            % rotate and shift vertices
            V = repmat(V,N_check,1) ;
            V_rot = R_blk*V + repmat(Z(:),1,size(V,2)) ;

            % each column of V_rot is the position of a vertex at each time in
            % T, where every three rows correspond to a single time; reshape
            % these to get the trajectories of each vertex
            V_trajs = nan(3*N_V,N_check) ;
            v_idx = 1 ;
            for idx = 1:N_V
                V_trajs(v_idx:v_idx+2,:) = reshape(V_rot(:,idx),3,[]) ;
                v_idx = v_idx + 3 ;
            end
            %V_trajs = reshape(V_rot(:),[],N_check) ;

            % interpolate vertices linearly (this is a coarse approximation but
            % hey it'll work for now)
            N_T_check = ceil((T(end) - T(1))/W.collision_check_time_discretization) ;
            T_check = linspace(T(1),T(end),N_T_check) ;
            if isempty(T_check)
                V_check = V_trajs ;
            else
                V_check = match_trajectories(T_check,T,V_trajs,'previous') ;
            end
            
            % V_check is now 3*N_V-by-N_T_check, where N_V is the number of
            % vertices; we pass the columns to the collision check loop:
        else
            V_check = [] ;
        end
    end
    
    function out = collision_check_loop(W,V_check)
        % initialize check to false optimistically
        out = false ;

        % V_check is now 3*N_V-by-N_T_check, where N_V is the number of
        % vertices; each column is a time slice, and every three rows is a
        % new vertex, so we need to reshape it to be 3-by-N
        V_check = reshape(V_check,3,[]) ;

        % FOR DEBUGGING
        % plot3(V_check(1,:),V_check(2,:),V_check(3,:),'r.')

        % iterate through obstacles and check for collisions
        obs_idx = 1 ;
        while ~out && obs_idx <= W.N_obstacles
            O_idx = W.obstacles{obs_idx} ;
            out = any(inpolyhedron(O_idx.collision_check_input_struct,V_check')) ;
            obs_idx = obs_idx + 1 ;
        end
    end

    %% goal check
    function out = goal_check(W,agent_info)
        z = agent_info.position ;
        switch W.goal_type
            case 'sphere'
                out = any(vecnorm(W.goal - z) <= W.goal_radius) ;
            case 'box'
                out = z(1,end) >= (W.bounds(2) - W.goal_radius) ;
            otherwise
                error('Set the world goal_type to ''sphere'' or ''box''')
        end
    end

    %% plotting
    function plot(W)
        hold_check = hold_switch() ;
        
        W.plot_start_goal_and_bounds()

        % plot obstacles (each obstacle takes care of its own plotting)
        for idx = 1:W.N_obstacles
            plot(W.obstacles{idx}) ;
        end
        
        hold_switch(hold_check) ;
    end

    function plot_at_time(W,t)
        hold_check = hold_switch() ;
        
        W.plot_start_goal_and_bounds() ;

        % plot obstacles (each obstacle takes care of its own plotting)
        for idx = 1:W.N_obstacles
            plot_at_time(W.obstacles{idx},t) ;
        end
        
        hold_switch(hold_check) ;
    end
    
    function plot_start_goal_and_bounds(W)
        if ~any(isinf(W.bounds))
            xlim([W.bounds(1), W.bounds(2)]);
            ylim([W.bounds(3), W.bounds(4)]);
            zlim([W.bounds(5), W.bounds(6)]);
        end

        if ~check_if_plot_is_available(W,'goal_zone')
            W.plot_goal() ;
        end
        
        if ~check_if_plot_is_available(W,'start_zone')
            W.plot_start()
        end
    end
    
    function plot_start(W)
        % plot start location
        if ~check_if_plot_is_available(W,'start')
            W.plot_data.start = plot3(W.start(1), W.start(2), W.start(3), 'ko') ;
        end
        
        % plot start zone
        if ~check_if_plot_is_available(W,'start_zone') && W.plot_start_zone_flag
            W.make_start_zone() ;
            S = W.start_zone ;
            W.plot_data.start_zone = patch('faces',S.faces,'vertices',S.vertices,...
                    'FaceColor',S.color,'EdgeColor',S.color,...
                    'FaceAlpha',S.opacity,'EdgeAlpha',S.opacity) ;
        end
    end
    
    function plot_goal(W)
        if isempty(W.goal_zone)
            W.make_goal_zone() ;
        end
        G = W.goal_zone ;
        
        switch W.goal_type
            case 'sphere'
                % plot goal location
                if ~check_if_plot_is_available(W,'goal')
                    W.plot_data.goal = plot3(W.goal(1), W.goal(2), W.goal(3), 'kx','LineWidth',1.5) ;
                end
                
                % plot ring around goal location
                xcirc = cos(linspace(0,2*pi)) ;
                ycirc = sin(linspace(0,2*pi)) ;
                r = W.goal_radius ;
                
                W.plot_data.goal_radius = plot3(W.goal(1)*ones(length(ycirc)),...
                    r*xcirc + W.goal(2),...
                    r*ycirc + W.goal(3),...
                    'k:','LineWidth',1.5) ;
                
                % plot goal zone as a sphere
                W.plot_data.goal_zone = surf(G.X,G.Y,G.Z,'FaceColor',G.color,'EdgeColor',G.color,...
                    'FaceAlpha',G.opacity,'EdgeAlpha',G.opacity) ;
            case 'box'
                W.plot_data.goal_box = patch('faces',G.faces,'vertices',G.vertices,...
                    'FaceColor',G.color,'EdgeColor',G.color,...
                    'FaceAlpha',G.opacity,'EdgeAlpha',G.opacity) ;
                goal_faces = zonotope_obstacle(W.bounds(2) - W.goal_zone.vertices(1,1),...
                                            W.bounds(4) - W.bounds(3),...
                                           W.bounds(6)-W.bounds(5),...
                                            [(W.goal_zone.vertices(1,1) + W.bounds(2))/2;
                                            (W.bounds(4) + W.bounds(3))/2;
                                            (W.bounds(5) + W.bounds(6))/2 ]);
                goal_faces.plot_face_opacity = G.opacity ;
                goal_faces.plot_face_color = G.color ;
                plot(goal_faces)
            otherwise
                error('Set the world goal_type to ''sphere'' or ''box''')
        end
    end
    
    function make_goal_zone(W)
        switch W.goal_type
            case 'sphere'
                [gx,gy,gz] = sphere(20) ;
                gx = W.goal_radius.*gx + W.goal(1) ;
                gy = W.goal_radius.*gy + W.goal(2) ;
                gz = W.goal_radius.*gz + W.goal(3) ;
                W.goal_zone.X = gx ;
                W.goal_zone.Y = gy ;
                W.goal_zone.Z = gz ;
            case 'box'
                B = W.bounds ;
                
                faces = [1 2 3 4] ;
                
                y_verts = [B(3) B(4) B(4) B(3)] ;
                z_verts = [B(5) B(5) B(6) B(6)] ;
                x_verts = (B(2) - W.goal_radius).*ones(1,4) ;
                vertices = [x_verts ; y_verts ; z_verts]' ;
                
                W.goal_zone.faces = faces ;
                W.goal_zone.vertices = vertices ;
            otherwise
                error('Set the world goal_type to ''sphere'' or ''box''')
        end
        
        W.goal_zone.color = W.plot_goal_zone_color ;
        W.goal_zone.opacity = W.plot_goal_zone_opacity ;
    end
    
    function make_start_zone(W)
        B = W.bounds ;
                
        faces = [1 2 3 4] ;
        
        y_verts = [B(3) B(4) B(4) B(3)] ;
        z_verts = [B(5) B(5) B(6) B(6)] ;
        x_verts = (B(1) + W.buffer_start).*ones(1,4) ;
        vertices = [x_verts ; y_verts ; z_verts]' ;
        
        W.start_zone.faces = faces ;
        W.start_zone.vertices = vertices ;
        W.start_zone.color = W.plot_start_zone_color ;
        W.start_zone.opacity = W.plot_start_zone_opacity ;
    end

    function set_obstacle_opacity(W,in,idxs)
        % W.set_obstacle_opacity(W,opacity_in,obs_indices)
        %
        % Pass in an opacity (alpha value) as either a scalar or a
        % 2-element vector. The face alphas are set using the first value
        % of the opacity input, and the edges with the second. If no
        % obstacle indices are passed in, the values are applied to all
        % obstacles.
        if nargin < 2
            in = 0.5 ;
        end
        if nargin < 3
            idxs = 1:W.N_obstacles ;
        end
        for idx = idxs
            W.obstacles{idx}.plot_face_opacity = in(1) ;
            if length(in) > 1
                W.obstacles{idx}.plot_edge_opacity = in(2) ;
            end
        end
    end

    function set_obstacle_color(W,c,idxs)
        % W.set_obstacle_color(color,obs_indices)
        %
        % Pass in a color as either a 1-by-3 (which changes the face
        % colors) or a 2-by-3 which changes the face and edge colors of the
        % obstacles specified by the obs_indices input. The color changes
        % are applied to all obstacles if no indices are provided.
        
        if nargin < 2
            c = [1 0 0] ;
        end
        if nargin < 3
            idxs = 1:W.N_obstacles ;
        end
        for idx = idxs
            W.obstacles{idx}.plot_face_color = c(1,:) ;
            if size(c,1) > 1
                W.obstacles{idx}.plot_face_color = c(2,:) ;
            end
        end
    end
end
end