classdef highway_cruising_agent < RTD_agent_2D
    % Class: highway_cruising_agent < RTD_agent_2D < agent
    %
    % Highway car closed loop dynamics. uses ode45 as integrator for
    % realism, have color plotting for wheel angle and body heading.
    %
    % based on real car model
    % states: x y phi vx delta
    
    properties
        lane
        lane_des
        % state limits
        %         % state indices
        speed_index = 4 ;
        
        % input limits 
        max_del_dot = 2; % rad/s
        max_accel = 4.0 ; % m/s^2
        
        wheel_plot_data={};
        wheel_color  = [130 130 130]/255;
        
        % integrator type, to allow for fixed time step integration
        integrator_type = 'ode45' ; % choose 'ode45' or 'ode4' or 'ode113'
        integrator_time_discretization = 0.01 ; % for ode4
        desired_initial_condition=[-500; 6; 0; 22; 0]; %start in the middle lane, speed will be changed by reset code
        
        footprint_vertices_for_plotting = [-2.4,-1.5,-1.5 0 0.3     2    2   2.4 2.4 2  2   0.3 0 -1.5 -1.5 -2.4 -2.4;
                -0.5,-0.5 -1   -1 -0.5    -0.3 -1   -1  1  1 0.3  0.5 1 1   0.5 0.5 -0.5];

        
    end
    
    methods
        %% constructor
        function A = highway_cruising_agent(varargin)
            % set up default superclass values
            name = 'highway_cruiser' ;
            
            default_footprint = [4.8 2];
            n_states = 5 ;
            n_inputs = 2 ; % two reference actually
            stopping_time = 50 ; % conservative estimate
            sensor_radius = 30 ;
            %             LLC = turtlebot_PD_LLC ; % doesnt need a controller
            
            % create agent
            A@RTD_agent_2D('name',name,...
                'footprint',default_footprint,...
                'n_states',n_states,'n_inputs',n_inputs,...
                'stopping_time',stopping_time,'sensor_radius',sensor_radius,varargin{:}) ;
            
            %draw F1 without wheels
           
            %             S.A.reset();
            
        end
        
        %% get agent info
        function agent_info = get_agent_info(A)
            % call superclass method
            agent_info = get_agent_info@RTD_agent_2D(A) ;
            agent_info.lane = A.lane;
            agent_info.lane_des = A.lane_des;
            
        end
        function plot(A)
            plot@RTD_agent_2D(A);
            A.plot_wheel_at_time(A.time(end))
        end
        function plot_at_time(A,t)
            plot_at_time@RTD_agent_2D(A,t);
            A.plot_wheel_at_time(t);
            z_t = match_trajectories(t,A.time,A.state) ;
            xlim([z_t(1)-20,z_t(1)+50]); ylim([-0.7, 12.7]);
        end
        function plot_wheel_at_time(A,t)
            wheel_size = [0.7 0.4];
            wheel = make_box(wheel_size);
            %rare , front
            wheel_position = [-2   -2  1.5 1.5
                -0.75 0.75 -0.75 0.75];
            wheel_vertices = [];
            for i = 1:4
                wheel_vertices = [wheel_vertices wheel+repmat(wheel_position(:,i),[1,5]) [NaN;NaN]];
            end
            % compute footprint for plot
            z_t = match_trajectories(t,A.time,A.state) ;
            p_t = z_t(A.position_indices) ;
            h_t = z_t(A.heading_index) ;
            delta_t = z_t(5);
            R_r = rotation_matrix_2D(h_t); %+ rotation_matrix_2D(delta_t)  ;
            V_ft = R_r*A.footprint_vertices_for_plotting + repmat(p_t,1,size(A.footprint_vertices_for_plotting,2));
            R_f = rotation_matrix_2D(h_t+delta_t);
            V_all = R_r*wheel_vertices + repmat(p_t,1,size(wheel_vertices,2)) ;
            %             V_front =
            %             V_fp = [V_rare V_front];
            
            for i = 1:4
                if i == 3 || i == 4
                    wheel_vert = V_all(:,6*i-5:6*i-1);
                    wheel_center = repmat( 0.5*(max(wheel_vert,[],2)+min(wheel_vert,[],2)),[1,5]);
                    origion_vert = wheel_vert - wheel_center;
                    V_all(:,6*i-5:6*i-1) = R_f * origion_vert + wheel_center;
                end
                %             fill(V_all(1,6*i-5:6*i-1),V_all(2,6*i-5:6*i-1),[80 80 80]/255)
                
            end
            if check_if_plot_is_available(A,'wheel_plot_data')
                for i = 1:4
                    A.plot_data.wheel_plot_data{i}.XData =  V_all(1,6*i-5:6*i-1) ;
                    A.plot_data.wheel_plot_data{i}.YData =  V_all(2,6*i-5:6*i-1);
                    uistack(A.plot_data.wheel_plot_data{i}, 'top')
                end
                
            else
                for i =1:4
                    h = fill( V_all(1,6*i-5:6*i-1), V_all(2,6*i-5:6*i-1),A.wheel_color) ;
                    A.plot_data.wheel_plot_data{i} = h ;
                    h.FaceAlpha = A.plot_footprint_opacity;
                    h.EdgeAlpha = A.plot_footprint_edge_opacity;
                end
            end
            if check_if_plot_is_available(A,'pretty_footprint')
                A.plot_data.pretty_footprint.Vertices = V_ft' ;
                uistack(A.plot_data.pretty_footprint, 'top')
            else
                % plot footprint
                fp_data = patch(V_ft(1,:),V_ft(2,:),A.plot_footprint_color,...
                    'EdgeColor',A.plot_footprint_edge_color,...
                    'FaceAlpha',A.plot_footprint_opacity,...
                    'EdgeAlpha',A.plot_footprint_edge_opacity) ;
                
            
                % save plot data
                A.plot_data.pretty_footprint = fp_data ;
            end
            
         
        end
        function reset(A,state)
            if nargin < 2
                start_lane = 1; %start on second lane
                y_ini = 4*start_lane + 2;
                A.lane = start_lane;
                A.lane_des = start_lane;
                A.desired_time = zeros(1,0);
                A.desired_input = zeros(2,0);
                A.desired_trajectory =zeros(2,0);
                start_v = 1;%always start with 1m/s
                if isprop(A,'desired_initial_condition')
                    A.desired_initial_condition(2) = y_ini;
                    A.desired_initial_condition(4) = start_v;
                    reset@RTD_agent_2D(A,[A.desired_initial_condition]) ;
                else
                    reset@RTD_agent_2D(A) ;
                end
            else
                reset@RTD_agent_2D(A,state) ;
            end
        end
        
        %% dynamics
        function dzdt = dynamics(A,t,z,T,U,Z)
            % handle no desired trajectory input
            if nargin < 6
                Z = [] ;
            end
            
            % extract the states
            y = z(2);
            psi =  z(A.heading_index);%3
            vx  =  z(A.speed_index);%4
            delta = z(5);
            
            
            vx_nom =interp1(T,U(1,:),t,'linear') ;
            y_nom = interp1(T,U(2,:),t,'linear') ; % speed
            if vx_nom < 0
                vx_nom = 0;
            end            
            
            Tc = 2;
            m = 1558;
            lf = 1.462884;
            lr = 1.405516;
            l = lf + lr;
            
            Cf = 1.432394487827058e+05;
            Cr = 2.214094963126969e+05;
            g = 9.80655;
            % get feedback control inputs
            
            ua = 5*(vx_nom-vx);
            delta_cmd =1*(y_nom-y)-0.1*psi-0.1*delta;
            delta_dot = 5*(delta_cmd-delta);
            
            
            kus = (m*g*lr/(l*Cf)-m*g*lf/(l*Cr));
            
            yr = tan(delta)*vx/(l+kus*vx^2/g);
            
            vy = yr*(lr-m*vx^2*lr/(Cr*l));
            
            
            % saturate the inputs
            if abs(delta_dot) > A.max_del_dot
                warning('delta saturated')
            end
            if abs(ua) > A.max_accel
                warning('throttle/brake saturated')
            end
            delta_dot = bound_values(delta_dot,A.max_del_dot) ;
            ua = bound_values(ua,A.max_accel) ;
            
            % calculate the derivatives
            dzdt = [vx*cos(psi)-vy*sin(psi);...
                vx*sin(psi)+vy*cos(psi);...
                yr;...
                -Tc/m*vx+Tc*ua;...
                delta_dot];
        end
        
        %% integrator options
        function [tout,zout] = integrator(A,fun,tspan,z0)
            switch A.integrator_type
                case 'ode45'
                    [tout,zout] = ode45(@(t,z) fun(t,z),tspan,z0(:)) ;
                case 'ode113'
                    [tout,zout] = ode113(@(t,z) fun(t,z),tspan,z0(:)) ;
                case {'ode4','RK4'}
                    dt = A.integrator_time_discretization ;
                    tout = tspan(1):dt:tspan(end) ;
                    if tout(end) ~= tspan(end)
                        tout = [tout, tspan(end)] ;
                    end
                    zout = ode4(@(t,z) fun(t,z),tout,z0(:)) ;
                otherwise
                    error('Please set A.integrator_type to either ode45 or ode4')
            end
            tout = tout(:)' ;
            zout = zout' ;
        end
    end
end