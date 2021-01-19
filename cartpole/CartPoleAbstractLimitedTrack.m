classdef (Abstract) CartPoleAbstractLimitedTrack < rl.env.MATLABEnvironment
    %RLENVCartPoleAbstractLimitedTrack: Abstract class for cart-pole inverted pendulum 
    %environments in MATLAB.
    
    % Copyright 2017-2018 The MathWorks Inc.    
    
    
    % Note: Do we want this to be more abstract for swing up
    % environments and rewards?
        
    properties
        % Acceleration due to gravity in m/s^2
        Gravity = 9.8
        
        % Mass of the cart
        MassCart = 1.0
        
        % Mass of the pole
        MassPole = 0.1
        
        % Half the length of the pole
        Length = 0.5
        
        % Max Force the input can appy
%         MaxForce = 10
        MaxForce = 1   
               
        % Sample time
        Ts = 0.02
        
        % Parameterization time
        Tf = 1
        
        % Angle at which to fail the episode
        ThetaThresholdRadians = 12 * pi/180 % 12 * pi/180
        
        % Distance at which to fail the episode
        XThreshold = 2.4    %2.4
        
        % The parameter for the controller
        k1 = 50;
        k2 = 5;
        
        % Define the error bound
        error_bound_min = [];
        error_bound_max = [];
        
        % define the sliced resolution and range
        % set the number of points we want to sample
        N_x = 5;
        N_x_dot = 10;
        N_theta = 10;
        N_theta_dot =10;
        % set the range of sampled space
        x_start = -2; x_end = 2;
        x_dot_start = -5; x_dot_end = 5;
        theta_start = 0; theta_end = 2*pi;
        theta_dot_start = -5; theta_dot_end = 5;

        
    end
    
    properties
        % system state [x,dx,theta,dtheta]'
        State = zeros(4,1)
    end
    
    properties(Access = protected)
        % Internal flag to store stale env that is finished
        IsDone = false
    end
    properties (Transient,Access = private)
        Visualizer = []
    end
    methods (Abstract,Access = protected)
        Force = getForce(this,action)
        updateActionInfo(this)
        Reward = getReward(this,x,force)
    end    
    methods                
        function this = CartPoleAbstractLimitedTrack(ActionInfo)
            ObservationInfo = rlNumericSpec([5 1]);
            ObservationInfo.Name = 'CartPole States';
            ObservationInfo.Description = 'sin(T0), cos(T0), Td0, X0, Xd0';
%             ObservationInfo.Description = 'x, theta, dx, dtheta';
            this = this@rl.env.MATLABEnvironment(ObservationInfo,ActionInfo); 
        end
        function set.State(this,state)
            validateattributes(state,{'numeric'},{'finite','real','vector','numel',4},'','State');
            this.State = double(state(:));
            notifyEnvUpdated(this);
        end
        function set.Length(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','Length');
            this.Length = val;
            notifyEnvUpdated(this);
        end
        function set.Gravity(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','Gravity');
            this.Gravity = val;
        end
        function set.MassCart(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','MassCart');
            this.MassCart = val;
        end
        function set.MassPole(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','MassPole');
            this.MassPole = val;
        end
        function set.MaxForce(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','MaxForce');
            this.MaxForce = val;
            updateActionInfo(this);
        end
        function set.Ts(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','Ts');
            this.Ts = val;
        end
        function set.ThetaThresholdRadians(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','ThetaThresholdRadians');
            this.ThetaThresholdRadians = val;
            notifyEnvUpdated(this);
        end
        function set.XThreshold(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','XThreshold');
            this.XThreshold = val;
            notifyEnvUpdated(this);
        end
%         function set.RewardForNotFalling(this,val)
%             validateattributes(val,{'numeric'},{'real','finite','scalar'},'','RewardForNotFalling');
%             this.RewardForNotFalling = val;
%         end
%         function set.PenaltyForFalling(this,val)
%             validateattributes(val,{'numeric'},{'real','finite','scalar'},'','PenaltyForFalling');
%             this.PenaltyForFalling = val;
%         end
        
        % we use ODE 45 to solve the step
        function [observation,reward,isdone,loggedSignals] = step(this,action)
            
            loggedSignals = [];
            
            % Get action
            force = getForce(this,action);            
            
            % Unpack state vector
            state = this.State;
            
%             % x = state(1);
%             x_dot = state(2);
%             theta = state(3);
%             theta_dot = state(4);
            
%             % Apply motion equations
%             costheta = cos(theta);
%             sintheta = sin(theta);
%             
%             totalmass = this.MassCart + this.MassPole;
%             polemasslength = this.MassPole*this.Length;
%             
%             temp = (force + polemasslength * theta_dot * theta_dot * sintheta) / totalmass;
%             thetaacc = (this.Gravity * sintheta - costheta* temp) / (this.Length * (4.0/3.0 - this.MassPole * costheta * costheta / totalmass));
%             xacc  = temp - polemasslength * thetaacc * costheta / totalmass;
%             
%             % Euler integration
%             observation = state + this.Ts.*[x_dot;xacc;theta_dot;thetaacc] + 0.5 * this.Ts.* this.Ts.*[xacc;0;thetaacc;0];

            Tspan = [0,this.Ts]; % define the timespan
            [t, state] = ode45(@(t,y) cartpole_dynamic(t,y,force,this), Tspan, state); % using ODE45 to solve the evolution of the system
%             cartpole_plot(this,t,state);
            observation = state(end,:).';
            this.State = observation;
            
            x = observation(1);
%             theta = observation(3);
            
%             isdone = abs(x) > this.XThreshold || abs(theta) > this.ThetaThresholdRadians;
            isdone = abs(x) > this.XThreshold;  % we only need to check the XTHreshold, theta can be anything
            this.IsDone = isdone;
            
            % Get reward
            reward = getReward(this,x,force);   
        end  
        
        function state_dot = cartpole_dynamic(t,state_in,force,this)
            % t is the input time
            % state_in is the input the state vector
            % force is the input the forcing, in th cartpole case, it's the force
            % exerted on the cart and we assume it's constant here
            % this is set as the input to give the geometory of the cart pole system

            % state_in contains four components, it's a 4x1 vector, state_in(1) is x, state_in(2) is theta,
            % state_in(3) is x_dot, state_in(4) is theta_dot

            % state_dot is the output,it's the derivative of state, and state_dot=f(state,F) is a
            % nonlinear function

            % obtain the geometry of the cart pole
            totalmass = this.MassCart + this.MassPole;
            polemasslength = this.MassPole*this.Length;

            % obtain the state
            theta = state_in(3);
            % x_dot = state_in(3);
            theta_dot = state_in(4);

            sintheta = sin(theta);
            costheta = cos(theta);


            temp = (force + polemasslength * theta_dot * theta_dot * sintheta) / totalmass;
            theta_acc = (this.Gravity * sintheta - costheta* temp) / (this.Length * (4.0/3.0 - this.MassPole * costheta * costheta / totalmass));
            x_acc  = temp - polemasslength * theta_acc * costheta / totalmass;

            % assign the value to y_dot
            state_dot = zeros(4,1);
%             state_dot(1) = state_in(3); 
%             state_dot(2) = state_in(4);
%             state_dot(3) = x_acc;
%             state_dot(4) = theta_acc;
            state_dot(1) = state_in(2); 
            state_dot(3) = state_in(4);
            state_dot(2) = x_acc;
            state_dot(4) = theta_acc;     
        end
        function [H1] = cartpole_plot(this,t,y)
            % the function is used to plot the four components of the state vector
            H1 = figure(1);
            subplot(2,2,1);
            plot(t,y(:,1));
            title('x');

            subplot(2,2,2);
            plot(t,y(:,2)/pi*180);
            title('theta');

            subplot(2,2,3);
            plot(t,y(:,3));
            title('x\_dot');

            subplot(2,2,4);
            plot(t,y(:,4));
            title('theta\_dot');

        end
        function initialState = reset(this)
            % Randomize the initial pendulum angle between (+- .05 rad)
            
            % Theta (+- .05 rad)
%             T0 = 2*0.05*rand - 0.05;  
            % in order to swing up, we put the initial condition as pi
            T0 = pi; %pi 
            % Thetadot
            Td0 = 0;
            % X 
            X0 = 0;
            % Xdot
            Xd0 = 0;
%             initialState= [T0;Td0;X0;Xd0];
%             initialState= [X0;T0;Xd0;Td0];
            initialState= [X0;Xd0;T0;Td0];
            this.State = initialState;
        end        
        function varargout = plot(this)
            % Visualizes the environment
            if isempty(this.Visualizer) || ~isvalid(this.Visualizer)
                this.Visualizer = CartPoleVisualizerLimitedTrackNew(this);
            else
                bringToFront(this.Visualizer);
            end
            if nargout
                varargout{1} = this.Visualizer;
            end
        end
    end
end
