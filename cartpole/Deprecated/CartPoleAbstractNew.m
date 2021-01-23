classdef (Abstract) CartPoleAbstractNew < MATLABEnvironment
    %RLENVCARTPOLEABSTRACT: Abstract class for cart-pole inverted pendulum 
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
        MaxForce = 10
               
        % Sample time
        Ts = 0.02
        
        % Angle at which to fail the episode
        ThetaThresholdRadians = 12 * pi/180
        
        % Distance at which to fail the episode
        XThreshold = 2.4
        
        % Reward each time step the cart-pole is balanced
        RewardForNotFalling = 1
        
        % Penalty when the cart-pole fails to balance
        PenaltyForFalling = -5 
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
        function this = CartPoleAbstractNew(ActionInfo)
            ObservationInfo = rlNumericSpec([4 1]);
            ObservationInfo.Name = 'CartPole States';
            ObservationInfo.Description = 'x, dx, theta, dtheta';
            this = this@MATLABEnvironment(ObservationInfo,ActionInfo); 
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
        function set.RewardForNotFalling(this,val)
            validateattributes(val,{'numeric'},{'real','finite','scalar'},'','RewardForNotFalling');
            this.RewardForNotFalling = val;
        end
        function set.PenaltyForFalling(this,val)
            validateattributes(val,{'numeric'},{'real','finite','scalar'},'','PenaltyForFalling');
            this.PenaltyForFalling = val;
        end
        function [observation,reward,isdone,loggedSignals] = step(this,action)
            
            loggedSignals = [];
            
            % Get action
            force = getForce(this,action);            
            
            % Unpack state vector
            state = this.State;
            % x = state(1);
            x_dot = state(2);
            theta = state(3);
            theta_dot = state(4);
            
            % Apply motion equations
            costheta = cos(theta);
            sintheta = sin(theta);
            
            totalmass = this.MassCart + this.MassPole;
            polemasslength = this.MassPole*this.Length;
            
            temp = (force + polemasslength * theta_dot * theta_dot * sintheta) / totalmass;
            thetaacc = (this.Gravity * sintheta - costheta* temp) / (this.Length * (4.0/3.0 - this.MassPole * costheta * costheta / totalmass));
            xacc  = temp - polemasslength * thetaacc * costheta / totalmass;
            
            % Euler integration
            observation = state + this.Ts.*[x_dot;xacc;theta_dot;thetaacc];

            this.State = observation;
            
            x = observation(1);
            theta = observation(3);
            
            isdone = abs(x) > this.XThreshold || abs(theta) > this.ThetaThresholdRadians;
            this.IsDone = isdone;
            
            % Get reward
            reward = getReward(this,x,force);   
        end  
        function initialState = reset(this)
            % Randomize the initial pendulum angle between (+- .05 rad)
            
            % Theta (+- .05 rad)
            T0 = 2*0.05*rand - 0.05;  
            % Thetadot
            Td0 = 0;
            % X 
            X0 = 0;
            % Xdot
            Xd0 = 0;
            
            initialState= [T0;Td0;X0;Xd0];
            this.State = initialState;
        end        
        function varargout = plot(this)
            % Visualizes the environment
            if isempty(this.Visualizer) || ~isvalid(this.Visualizer)
                this.Visualizer = rl.env.viz.CartPoleVisualizer(this);
            else
                bringToFront(this.Visualizer);
            end
            if nargout
                varargout{1} = this.Visualizer;
            end
        end
    end
end