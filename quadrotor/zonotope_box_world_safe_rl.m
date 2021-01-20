classdef zonotope_box_world_safe_rl < zonotope_box_world
    %% properties
    properties
        % number of points to use to create obstacle observations (this
        % value is cubed to become a 3-D grid of points)
        N_obstacle_observation_points = 2 ;
        
        % sphere of obstacle observation points
        d_obstacle_observation_points = 5 ; % radius of this sphere
        obstacle_observation_points = [] ; % store the points in this property
        ob_dirs=[];
        el_mat
        az_mat
    end
    
    %% methods
    methods
        %% constructor
        function W = zonotope_box_world_safe_rl(varargin)
            W@zonotope_box_world(varargin{:}) ;
            
            % create obstacle observation points
            d_max = W.d_obstacle_observation_points ;
            N_p = W.N_obstacle_observation_points ;
            p_vec = linspace(-d_max, +d_max, N_p) ;
            [P1,P2,P3] = ndgrid(p_vec,p_vec,p_vec) ;
            P = [P1(:), P2(:), P3(:)]' ;
            
            % make these points into a sphere
            d = vecnorm(P) ;
            d_log = d <= d_max ;
            P_ob = P(:,d_log) ;
            
            % store the points
            W.obstacle_observation_points = P_ob ;
            
            %sample az angle and el angle. make observation pts, see paper appendix for details:
            % https://arxiv.org/abs/2011.08421
            az = linspace(-pi/4,pi/4,9);
            az = az(2:end-1); % start and end same
            az = [-pi/2 az pi/2];
            el = linspace(-pi/3,pi/3,9);
            el = el(4:end-3); %top and bottom shouldnt be repeated, just remove

            az_mat = repmat(az,[length(el),1]);
            az_mat = az_mat(:);
            az_mat = [az_mat; 0; 0];
            el_mat = repmat(el',[1,length(az)]);
            el_mat = el_mat(:);
            el_mat = [el_mat; pi/2;-pi/2];
            [x,y,z] = sph2cart(az_mat,el_mat,ones(length(el_mat),1));

            W.ob_dirs = [x';y';z']; %unit vectors pointing at the direction of the desired place
            W.el_mat = el_mat;
            W.az_mat = az_mat;
            
        end
        
        %% get observation
        function ob = get_ob(W,agent_info)
            % ob = W.get_ob(agent_info)
            %
            % Create an Observation output for RL. For the quadrotor, we
            % return the following:
            %   1. the quadrotor's initial velocity
            %   2. the unit vector pointing towards the global goal
            %   3. a sephere of directions that shows the distance to
            %   obstacles, angles picked very carefully
            
            % create observation with initial velocity and acceleration
            % (which are part of the trajectory parameterization, so it's
            % useful for the RL agent to learn them)
            ob = [agent_info.velocity(:,end)] ;
            
            % create a unit vector pointing towards the global goal
            x_cur = agent_info.position(:,end) ;
            g = W.goal +[1.5;0;0];
            u = make_unit_length(g - x_cur) ;
            ob = [ob ; u] ;
            
            % 4::
            % prep arrays
            N_ob = size(W.ob_dirs,2) ;
            D = ones(1,N_ob)*agent_info.sensor_radius;
            for obser_idx = 1:N_ob
                obser_dir = W.ob_dirs(:,obser_idx);
                for obst_idx = 1:W.N_obstacles
                    o = W.obstacles{obst_idx} ;% including boundaries
                    [flag, dist] = rayBoxIntersection(x_cur, obser_dir, min(o.vertices,[],1), max(o.vertices,[],1));
                    if flag && dist > 0
                       D(obser_idx) = min(dist,D(obser_idx));
                    end
                end
            end
%               uncomment the following to visualize it during manual control/training/evaluation/
%             [x,y,z] = sph2cart(W.az_mat(:),W.el_mat(:),D(:));
%             scatter3(x+ones(N_ob,1)*x_cur(1),y+ones(N_ob,1)*x_cur(2),z+ones(N_ob,1)*x_cur(3)); axis equal
%             end
       
            ob = [ob ; D(:)] ;
        end
        
        %% get reward
        function total_rew = getRew(W,agent_info,observation)
            % reward is designed so it is mostly negative, getting to the goal gets you a big reward outside
            rew_v = 0.5*dot(agent_info.velocity(:,end),observation(4:6))-2; %project velocity in the direction of goal
            %0.5 to -4.5 since velocity max is 5;
            dist_vec = sort(observation(7:end));
            rew_obst = atan(mean(dist_vec(1:8+5)))^4-4;%range -4 ~ 0.6; % there are 8 observations facing the top and bottom wall, so to get flying information in front, take 5 closest points in front.
            %Make sure all 13 are not too close
            total_rew = rew_obst + rew_v; %+;% - agent_info.replace_distance/7;
            
        end
    end
end
