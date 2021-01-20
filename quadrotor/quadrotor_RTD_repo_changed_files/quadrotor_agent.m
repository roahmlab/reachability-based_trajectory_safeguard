classdef quadrotor_agent < rigid_body_agent_SE3
properties
    % the default properties are for an AscTec Hummingbird quadrotor
    %
    % http://www.asctec.de/downloads/public/HB_AscTec-Hummingbird_safety-data-sheet.pdf
    % https://journals.sagepub.com/doi/pdf/10.5772/59618
    % https://digitalrepository.unm.edu/cgi/viewcontent.cgi?referer=https://www.google.com/&httpsredir=1&article=1189&context=ece_etds
    %
    % The third link is to William Neeley's master's thesis. See Ch 8
    % starting on pg. 113 of the PDF for rotor speed constraints.

    % trajectory info
    T_des
    U_des
    Z_des
    
    % current state estimates
    a_est
    j_est
    
    % move method -- choose 'direct' or 'integrator'
    % the direct mode takes in a reference trajectory and "backs out" the
    % robot's required states along that trajectory using differential
    % flatness, whereas the integrator mode uses the robot's low-level
    % controller (the A.LLC property), the robot dynamics, and an ode
    % solver to simulate the robot's closed-loop high-fidelity model when
    % tracking a given trajectory plan; the default integrator is ode45
    move_method = 'integrator' ;
    
    % structures for physical properties
    body_dimensions = [0.5 0.5 0.05] ;
    body_distance_to_rotor = 0.54/2 ; % m from center of body_mass
    
    rotor_thrust_coefficient = 1.5e-7
    rotor_moment_coefficient = 3.75e-9 ;
    rotor_hub_displacements
    rotor_speed_squared_to_inputs_matrix % matrix used to map rotor speeds to inputs
    rotor_speed_squared_to_inputs_matrix_inverse
    rotor_min_speed_squared = 1100^2 ; % rpm
    rotor_max_speed_squared = 8600^2 ; % rpm
    
    % camera parameters
    camera_follow_distance = 3 ; % distance from camera to quadrotor
    camera_direction_offset = [-3.0; 0.0; 1.5] ;  % direction of camera rel. to QR
    camera_target_offset = [ 2.0; 0.0; 0.0] ;  % target of camera rel. to QR
    camera_view_style = 'none' ; % none, above, behind, onboard
    
    % plot params
    plot_body_face_color = [0.7 0.7 1] ;
    plot_body_edge_color = [0 0 1] ;
    plot_body_opacity = 1 ;
end
methods
%% constructor
    function A = quadrotor_agent(varargin)
        % agent = quadrotor_agent()
        %
        % Creates a simulator agent modeled after an AscTec Hummingbird
        
        % set up default mass and inertia
        m = 0.547 ; % kg
        J = diag([0.0033 0.0033 0.0058]) ;
        A@rigid_body_agent_SE3('body_mass',m,...
            'body_inertia_matrix',J,...
            'integrator_method','ode45',...
            varargin{:}) ;
        
        % set up default vectors from rotor centers to QR center of body_mass
        l = norm(A.body_dimensions(1:2)./2) ;
        h = A.body_dimensions(3)./2 ;
        A.rotor_hub_displacements = [l -l -l l ; l l -l -l ; h h h h] ;
                                           
        % set up default thrust/rotor speed conversion matrices
        b = A.rotor_thrust_coefficient ;
        k = A.rotor_moment_coefficient ;
        l = A.body_distance_to_rotor ;
        bl = b*l ;
        A.rotor_speed_squared_to_inputs_matrix = [b b b b ; 
                                          0 bl 0 -bl ;
                                          -bl 0 bl 0 ;
                                          k -k k -k] ;
        A.rotor_speed_squared_to_inputs_matrix_inverse = inv(A.rotor_speed_squared_to_inputs_matrix) ;
        
        % default sensor radius
        A.sensor_radius = max(A.sensor_radius,10) ; % m
        
        % set up collision check data
        A.make_collision_check_input_data() ;
        
        % set up plot data
        A.make_plot_input_data() ;
        A.plot_data.body = [] ;
        A.plot_data.desired_trajectory = [] ;
        A.plot_frame_scale = 0.07 ;
        A.plot_frame_linewidth = 3 ;
        A.animation_set_axes_flag = false ;
        
        % set up default low-level controller
        if isempty(A.LLC)
            A.LLC = Mellinger_LLC() ;
        end
        
        % fix number of states and inputs
        A.n_states = 9 ;
        A.n_inputs = 4 ;
        
        A.reset() ;
    end
    
%% reset
    function reset(A,position,attitude)
        if nargin < 2
            position = zeros(A.n_states,1) ;
        end
        
        if nargin < 3
            attitude = eye(3) ;
        end
        
        reset@rigid_body_agent_SE3(A,position,attitude) ;
        
        % reset desired trajectory
        A.T_des = [] ;
        A.U_des = [] ;
        A.Z_des = [] ;
        
        % reset state estimates
        A.a_est = zeros(3,1) ;
        A.j_est = zeros(3,1) ;
    end
    
%% get agent info
    function agent_info = get_agent_info(A)
        agent_info = get_agent_info@rigid_body_agent_SE3(A) ;
        agent_info.position = A.state(A.position_indices,:) ;
        agent_info.velocity = A.state(A.velocity_indices,:) ;
        agent_info.acceleration = A.a_est ;
        agent_info.attitude = A.attitude ;
        agent_info.collision_check_input_data.faces = A.collision_check_input_data.faces ;
        agent_info.collision_check_input_data.vertices = A.collision_check_input_data.vertices ;
        agent_info.body_dimensions = A.body_dimensions ;
    end
    
%% move and commit move data
    function move(A,t_move,T_ref,U_ref,Z_ref)
        % move(A,t_move,T_ref,U_ref,Z_ref)
        %
        % Unlike the superclas method, this saves the desired trajectory,
        % and lets you directly treat a provided (i.e., nominal or
        % reference) trajectory as the trajectory executed by the quadrotor
        if nargin < 5
            Z_ref = [] ;
        end
        
        switch A.move_method
            case 'direct'
                A.move_direct(t_move,T_ref,U_ref,Z_ref)
            case 'integrator'
                move@rigid_body_agent_SE3(A,t_move,T_ref,U_ref,Z_ref)
            otherwise
                error(['Please set the quadrotor_agent move_method property ',...
                    'to either direct or integrator'])
        end
    end
    
    function move_direct(A,t_move,T_ref,U_ref,Z_ref)
        if nargin < 5
            error(['The quadrotor_agent requires a reference trajectory ',...
                'Z_ref as an input to its move_direct method, since this ',...
                'reference is used to back out the full state of the robot.'])
        end
        
        [T_used,U_used,Z_used] = A.move_setup(t_move,T_ref,U_ref,Z_ref) ;
        
        % iterate through T_used to get the states Z
        N_t = length(T_used) ;
        z_out = nan(A.n_states,N_t) ;
        R_out = nan(3,3,N_t) ;
        for idx = 1:N_t
            t_idx = T_used(idx) ;
            [y,R] = get_quadrotor_state_using_diff_flatness(t_idx,T_used,Z_used) ;
            z_out(:,idx) = y ;
            R_out(:,:,idx) = R ;
        end
        
        % commit move data
        A.commit_move_data(T_used,z_out,R_out,T_used,U_used,Z_used) ;
    end
    
    function commit_move_data(A,tout,zout,Rout,T_used,U_used,Z_used)
        % commit move data as the regular rigid body
        commit_move_data@rigid_body_agent_SE3(A,tout,zout,Rout,T_used,U_used) ;
        
        % save the time
        if ~isempty(A.T_des)
            t_cur = A.T_des(end-1) ;
        else
            t_cur = 0 ;
        end
        A.T_des = [A.T_des, T_used(2:end) + t_cur, nan(1,1)] ;
        
        % save the input
        n_U = size(U_used,1) ;
        A.U_des = [A.U_des, U_used(:,2:end), nan(n_U,1)] ;
        
        % save the desired trajectory
        if ~isempty(Z_used)
            n_Z = size(Z_used,1) ;
            A.Z_des = [A.Z_des, Z_used(:,2:end), nan(n_Z,1)] ;        
        end
    end
    
%% stop
    function stop(A,~)
        % execute a braking maneuver
        
        % get current speed and acceleration
        v_0 = A.state(A.velocity_indices,end) ;
        a_0 = A.a_est - A.gravity_direction*A.gravity_acceleration ;
        
        % create spline
        v_peak = zeros(3,1) ;
        t_plan = 0.5 ;
        t_peak = 1 ;
        t_total = 3 ;
        dt = 0.05 ;
        t_extra = 2 ;
        [T,Z,~] = generate_spline_peak_speed(v_0,a_0,v_peak,t_plan,...
                                             t_peak,t_total,dt,t_extra) ;
                                         
        % get input
        T_in = T ;
        U_in = zeros(4,size(T,2)) ;
        Z_in = Z ;
        Z_in(1:3,:) = Z_in(1:3,:) + A.state(A.position_indices,end) ;
        
        % move agent
        A.move(t_peak,T_in,U_in,Z_in) ;
    end
    
%% dynamics
    function zd = dynamics(A,t,z,R,T_ref,U_ref,Z_ref)
        % call low-level controller
        U = A.LLC.get_control_inputs(A,t,z,T_ref,U_ref,Z_ref,R,A.a_est,A.j_est) ;
        
        % compute desired rotor speeds
        W = A.rotor_speed_squared_to_inputs_matrix_inverse*U ;
        
        % bound rotor speeds
        W = saturate(W,A.rotor_min_speed_squared,A.rotor_max_speed_squared) ;
        
        % compute new thrust and moments
        U_in = A.rotor_speed_squared_to_inputs_matrix*W ;
        
        % convert input to body force and moments (in global frame)
        F_in = U_in(1)*R*A.e3 ;
        M_in = U_in(2:4) ;
        
        % pass new forces and moments to rigid body dynamics as
        % zero-order-hold inputs for the integration time step duration;
        % note that the input T_ref_in is made longer than a single
        % integration time step to prevent NaNs from occurring if the rigid
        % body agent's dynamics need to interpolate T_ref_in and U_ref_in
        T_ref_in = [t, t + 2*A.integrator_time_discretization] ;
        U_ref_in = repmat([F_in;M_in],1,2) ;
        zd = dynamics@rigid_body_agent_SE3(A,t,z,R,T_ref_in,U_ref_in,Z_ref) ;
        
        % update acceleration and jerk estimates
        a_old = A.a_est ;
        A.a_est = F_in/A.body_mass + A.gravity_direction*A.gravity_acceleration ;
        A.j_est = (a_old - A.a_est)./A.integrator_time_discretization ;
    end
    
%% plotting
    function plot(A,~)
        plot@rigid_body_agent_SE3(A) ;
        
        hold_check = ishold ;
        
        if ~hold_check
            hold on ;
        end
        
%         % plot reference traj
%         if ~isempty(A.Z_des)
%             Z = A.Z_des(1:3,:) ;
%             if check_if_plot_is_available(A,'desired_trajectory')
%                 A.plot_data.desired_trajectory.XData = Z(1,:) ;
%                 A.plot_data.desired_trajectory.YData = Z(2,:) ;
%                 A.plot_data.desired_trajectory.ZData = Z(3,:) ;
%             else
%                 ref_traj_data = plot3(Z(1,:),Z(2,:),Z(3,:),'b--') ;
%                 A.plot_data.desired_trajectory = ref_traj_data ;
%             end
%         end
        
        % plot body
        R_c = A.attitude(:,:,end) ;
        z_c = A.state(:,end) ;
        A.plot_body(R_c,z_c) ;
        
        if ~hold_check
            hold off ;
        end
    end
    
    function plot_at_time(A,t)
        if nargin < 2
            t = 0 ;
        end

        % plot rigid body
        plot_at_time@rigid_body_agent_SE3(A,t) ;
        
        % plot quadrotor body
        [R_t,z_t] = A.get_state_and_attitude_at_time(t) ;
        A.plot_body(R_t,z_t) ; 
    end
    
    function plot_body(A,R,z)
        % get position
        p = z(A.position_indices) ;
        
        % get faces and vertices for plotting body
        FB1 = A.plot_input_data.body_faces.top_and_bottom ;
        FB2 = A.plot_input_data.body_faces.sides ;
        VB = A.plot_input_data.body_vertices ;
        
        FR = A.plot_input_data.rotor_faces ;
        VR = A.plot_input_data.rotor_vertices ;
        
        % move body and rotors to current position and attitude
        VB = (R*VB')' + repmat(p(:)',A.plot_input_data.body_N_vertices,1);
        VR = (R*VR')' + repmat(p(:)',A.plot_input_data.rotor_N_vertices,1);
        
        if check_if_plot_is_available(A,'body')
            A.plot_data.body.Vertices = VB ;
            A.plot_data.sides.Vertices = VB ;
            A.plot_data.rotor.Vertices = VR ;
        else
            body_patch_data =  patch('Faces',FB1,'Vertices',VB,...
                               'FaceAlpha',A.plot_body_opacity,...
                               'FaceColor',A.plot_body_face_color,...
                               'EdgeColor',A.plot_body_edge_color,...
                               'EdgeAlpha',A.plot_body_opacity) ;
            sides_patch_data = patch('Faces',FB2,'Vertices',VB,...
                               'FaceAlpha',A.plot_body_opacity,...
                               'FaceColor',A.plot_body_face_color,...
                               'EdgeColor',A.plot_body_edge_color,...
                               'EdgeAlpha',A.plot_body_opacity) ;
            rotor_patch_data = patch('Faces',FR,'Vertices',VR,...
                               'FaceAlpha',A.plot_body_opacity,...
                               'FaceColor',A.plot_body_face_color,...
                               'EdgeColor',A.plot_body_edge_color,...
                               'EdgeAlpha',A.plot_body_opacity) ;
            A.plot_data.body = body_patch_data ;
            A.plot_data.sides = sides_patch_data ;
            A.plot_data.rotor = rotor_patch_data ;
        end
        
        % set plot view
        A.set_plot_view(p) ;
    end
    
    function set_plot_view(A,p)
        switch A.camera_view_style
            case 'behind'
                cam_dir = A.camera_direction_offset ;
                cam_dir = cam_dir/norm(cam_dir) ;

                % go behind robot according to its current speed
                p_cam = p + A.camera_follow_distance*cam_dir ;

                % set camera
                ax = gca ;
                ax.CameraPosition = p_cam' ;
                ax.CameraTarget = p + A.camera_target_offset ;
                ax.Projection = 'perspective' ;

                ax_x = [-40 40] + p(1) ;
                ax_y = [-10 10] + p(2) ;
                ax_z = [-12 10] + p(3) ;
                ax_vals = [ax_x ax_y ax_z] ;

                axis(ax_vals(:)') ;
            case 'above'
                cam_dir = [0;0;1] ;

                % go above robot
                p_cam = p + A.camera_follow_distance*cam_dir ;

                % set camera
                ax = gca ;
                ax.CameraPosition = p_cam' ;
                ax.CameraTarget = p + [0;0;-A.camera_follow_distance] ;
                ax.Projection = 'perspective' ;

                ax_x = 8*[-1 1] + p(1) ;
                ax_y = 8*[-1 1] + p(2) ;
                ax_z = [-12 12] + p(3) ;
                ax_vals = [ax_x ax_y ax_z] ;

                axis(ax_vals(:)') ;
            case 'onboard'
                cam_dir = R*[1;0;0] ;

                % go on robot
                p_cam = p + [-0.2;0;0.1] ;

                % set camera
                ax = gca ;
                ax.CameraPosition = p_cam' ;
                ax.CameraTarget = p + 0.1*A.camera_follow_distance*cam_dir ;
                ax.Projection = 'perspective' ;

                ax_x = [-40 40] + p(1) ;
                ax_y = [-5 5] + p(2) ;
                ax_z = [-12 12] + p(3) ;
                ax_vals = [ax_x ax_y ax_z] ;

                axis(ax_vals(:)') ;
        end
    end
    
    function make_plot_input_data(A)
        % create quadrotor body
        [FB1,FB2,VB] = make_quadrotor_body() ;
        
        % get rotor centers
        cx = [0.17, 0, -0.17, 0] ;
        cy = [0, 0.17, 0, -0.17] ;
        
        % create quadrotor rotors
        [rx,ry,rz] = cylinder([0.1 0.1],20) ;
        rz(1,:) = -0.001 ;
        rz(2,:) = 0.005 ;
        % rz(3,:) = 0.003 ;
        
        FR = zeros(1,4) ;
        VR = [] ;
        
        for idx = 1:4
            % get coordinates of new rotor
            rx_idx = rx + cx(idx) ;
            ry_idx = ry + cy(idx) ;
            
            % convert new rotor to patch
            [f,v,~] = surf2patch(rx_idx,ry_idx,rz) ;
            
            % add vertices and faces
            VR = [VR;v] ;
            FR = [FR ; (f + max(FR(:)))] ;            
        end
        FR = FR(2:end,:) ;
        
        % add to agent object
        A.plot_input_data.rotor_faces = FR ;
        A.plot_input_data.rotor_vertices = VR ;
        A.plot_input_data.body_faces.top_and_bottom = FB1 ;
        A.plot_input_data.body_faces.sides = FB2 ;
        A.plot_input_data.body_vertices = VB ;
        A.plot_input_data.body_N_vertices = size(VB,1) ;
        A.plot_input_data.rotor_N_vertices = size(VR,1) ;
    end
    
    function make_collision_check_input_data(A)
        R = eul2rotm([0 0 pi/2],'XYZ') ;
        [F,V] = make_cuboid_for_patch(A.body_dimensions(1),...
                                      0.05*A.body_dimensions(2),...
                                      0.05*A.body_dimensions(3),...
                                      [0;0;-0.25*A.body_dimensions(3)]) ;
        V = [V ; (R*V')'] ;
        A.collision_check_input_data.faces = [F ; (F+8)] ;
        A.collision_check_input_data.vertices = V ;
    end
    
    %% utility
    function [v_stats, a_stats] = report_speed_and_accel_stats(A)
        % get speed stats
        v = vecnorm(A.state(A.velocity_indices,:)) ;
        v_avg = mean(v) ;
        v_min = min(v) ;
        v_max = max(v) ;
        
        v_stats = [v_min,v_avg,v_max] ;
        
        % numerically differentiate to get accel
        t = A.time ;
        a = diff(v,1,2)./diff(t) ;
        
        a_avg = mean(a,'omitnan') ;
        a_min = min(a,[],'omitnan') ;
        a_max = max(a,[],'omitnan') ;
        
        a_stats = [a_min,a_avg,a_max] ;
    end
end
end