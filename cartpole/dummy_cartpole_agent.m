classdef dummy_cartpole_agent < RTD_agent_2D
    % Class: dummy_cartpole_agent < RTD_agent_2D < agent
    %
    % Dummy agent is called by advanced_move
    %

    properties    
        Gravity = 9.8;
        MassCart = 1.0;
        MassPole = 0.1;
        Length = 0.5;
        MaxForce = 1;
        Ts = 0.1;
        Tf = 2;
        ThetaThresholdRadians = 12 * pi/180; % 12 * pi/180
        XThreshold = 3.5;    %2.4
        k1= 50;
        k2 =50;
        max_force = -inf;
        max_x_acc = -inf;

        
        plottor;
        plottor_flag = 0;
        integrator_type = 'ode45' ;
    end
    
    methods
        %% constructor
        function A = dummy_cartpole_agent(varargin)
            % set up default superclass values
            name = 'dummy_cartpole' ;
            
            default_footprint = 1;
            n_states = 4 ;
            n_inputs = 1 ; % two reference actually
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
            %this is rewritting the parent function instead of inheriting
            
%             agent_info = get_agent_info@RTD_agent_2D(A) ;
%             agent_info.lane = A.lane;
%             agent_info.lane_des = A.lane_des;
            
        end

        function reset(A,state)
            T0 = pi; %pi
            % Thetadot
            Td0 = 0;
            % X
            X0 = 0;
            % Xdot
            Xd0 = 0;
%             this.a = 0;
            
            if exist('state','var')
                X0 = state(1);
                Xd0 = state(2);
                T0 = state(3);
                Td0 = state(4);
%                 this.a = ic(5);
                state = [X0;Xd0;T0;Td0];
                reset@RTD_agent_2D(A,state) ;
            else
                state = [X0;Xd0;T0;Td0];
                reset@RTD_agent_2D(A,state) ;
            end
            

            
        end
        
        %% emergency stop
        function stop(A,t_stop)
            
            %             if nargin < 2
            %                 t_stop = A.stopping_time ;
            %             end
            %
            %             % get the current speed
            %             v = A.state(A.speed_index,end) ;
            %
            %             % check how long it will take to come to a stop and make the
            %             % stopping time vector
            %             t_req_to_stop = v/A.max_accel ;
            %             T_input = [0, max(t_req_to_stop,t_stop)] ;
            %
            %             % generate the input and desired trajectory
            %             U_input = zeros(2,2) ;
            %             Z_input = [repmat(A.state(1:3,end),1,2); zeros(1,2)] ;
            %
            %             % call move method to perform stop
            %             % A.LLC.gains = A.LLC.stop_gains ;
            %             A.move(t_stop,T_input,U_input,Z_input) ;
            %
            %             % reset the default gains after stopping
            %             % A.LLC.gains = A.LLC.default_gains ;
        end
        %% dynamics
        function state_dot = dynamics(A,t,z,T,U,Z)
            % handle no desired trajectory input
%             function state_dot = cartpole_dynamic_new(this,t,state_in,T,U,Z);
            this = A;
            state_in = z;
            % t is the input time
            % state_in is the input the state vector
            % force is the input the forcing, in th cartpole case, it's the force
            % exerted on the cart and we assume it's constant here
            % this is set as the input to give the geometory of the cart pole system
            
            % state_in contains four components, it's a 4x1 vector, state_in(1) is x, state_in(2) is theta,
            % state_in(3) is x_dot, state_in(4) is theta_dot
            
            % state_dot is the output,it's the derivative of state, and state_dot=f(state,F) is a
            % nonlinear function
            x_nom  = interp1(T,Z,t,'linear') ;
            vx_nom  = interp1(T,U,t,'linear') ;
            % obtain the geometry of the cart pole
            totalmass = this.MassCart + this.MassPole;
            polemasslength = this.MassPole*this.Length;
            
            % obtain the state
            theta = state_in(3);
            % x_dot = state_in(3);
            theta_dot = state_in(4);
            
            %llc
            x = state_in(1);
            x_dot = state_in(2);
            %             Force = k1*(action-x) - k2*x_dot;
            %LLC
            force = this.k1*(x_nom-x)+ this.k2*(vx_nom-x_dot);
%             if abs(force)> this.max_force
%                 this.max_force = abs(force);
%                 this.T = T;
%                 this.Z = Z;
%                 this.U = U;
%             end
            
            sintheta = sin(theta);
            costheta = cos(theta);
            
            
            temp = (force + polemasslength * theta_dot * theta_dot * sintheta) / totalmass;
            theta_acc = (this.Gravity * sintheta - costheta* temp) / (this.Length * (4.0/3.0 - this.MassPole * costheta * costheta / totalmass));
            x_acc  = temp - polemasslength * theta_acc * costheta / totalmass;
            if abs(x_acc)>this.max_x_acc
                this.max_x_acc = abs(x_acc);
            end
            % assign the value to y_dot
            state_dot = zeros(4,1);
            %             state_dot(1) = state_in(3);
            %             state_dot(2) = state_in(4);
            %             state_dot(3) = x_acc;
            %             state_dot(4) = theta_acc;
            state_dot(1) = state_in(2) - 0.2*  state_in(2);
            state_dot(3) = state_in(4) - 0.2 * state_in(4) ;
            state_dot(2) = x_acc;
            state_dot(4) = theta_acc;
%         end
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