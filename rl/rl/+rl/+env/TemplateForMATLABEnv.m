classdef TemplateForMATLABEnv < rl.env.MATLABEnvironment
    %TEMPLATEFORMATLABENV: Template for defining custom environment in MATLAB.    
    
    %% Properties (set properties' attributes accordingly)
    properties
        % Specify and initialize environment's necessary properties    
        % Acceleration due to gravity in m/s^2
        Gravity = 9.8
        
        % Mass of the cart
        CartMass = 1.0
        
        % Mass of the pole
        PoleMass = 0.1
        
        % Half the length of the pole
        HalfPoleLength = 0.5
        
        % Max Force the input can apply
        MaxForce = 10
               
        % Sample time
        Ts = 0.02
        
        % Angle at which to fail the episode (radians)
        AngleThreshold = 12 * pi/180
        
        % Distance at which to fail the episode
        DisplacementThreshold = 2.4
        
        % Reward each time step the cart-pole is balanced
        RewardForNotFalling = 1
        
        % Penalty when the cart-pole fails to balance
        PenaltyForFalling = -10 
    end
    
    properties
        % Initialize system state [x,dx,theta,dtheta]'
        State = zeros(4,1)
    end
    
    properties(Access = protected)
        % Initialize internal flag to indicate episode termination
        IsDone = false        
    end

    %% Necessary Methods
    methods              
        % Contructor method creates an instance of the environment
        % Change class name and constructor name accordingly
        function this = TemplateForMATLABEnv()
            % Initialize Observation settings
            ObservationInfo = rlNumericSpec([4 1]);
            ObservationInfo.Name = 'CartPole States';
            ObservationInfo.Description = 'x, dx, theta, dtheta';
            
            % Initialize Action settings   
            ActionInfo = rlFiniteSetSpec([-1 1]);
            ActionInfo.Name = 'CartPole Action';
            
            % The following line implements built-in functions of RL env
            this = this@rl.env.MATLABEnvironment(ObservationInfo,ActionInfo);
            
            % Initialize property values and pre-compute necessary values
            updateActionInfo(this);
        end
        
        % Apply system dynamics and simulates the environment with the 
        % given action for one step.
        function [Observation,Reward,IsDone,LoggedSignals] = step(this,Action)
            LoggedSignals = [];
            
            % Get action
            Force = getForce(this,Action);            
            
            % Unpack state vector
            XDot = this.State(2);
            Theta = this.State(3);
            ThetaDot = this.State(4);
            
            % Cache to avoid recomputation
            CosTheta = cos(Theta);
            SinTheta = sin(Theta);            
            SystemMass = this.CartMass + this.PoleMass;
            temp = (Force + this.PoleMass*this.HalfPoleLength * ThetaDot^2 * SinTheta) / SystemMass;

            % Apply motion equations            
            ThetaDotDot = (this.Gravity * SinTheta - CosTheta* temp) / (this.HalfPoleLength * (4.0/3.0 - this.PoleMass * CosTheta * CosTheta / SystemMass));
            XDotDot  = temp - this.PoleMass*this.HalfPoleLength * ThetaDotDot * CosTheta / SystemMass;
            
            % Euler integration
            Observation = this.State + this.Ts.*[XDot;XDotDot;ThetaDot;ThetaDotDot];

            % Update system states
            this.State = Observation;
            
            % Check terminal condition
            X = Observation(1);
            Theta = Observation(3);
            IsDone = abs(X) > this.DisplacementThreshold || abs(Theta) > this.AngleThreshold;
            this.IsDone = IsDone;
            
            % Get reward
            Reward = getReward(this);
            
            % (optional) use notifyEnvUpdated to signal that the 
            % environment has been updated (e.g. to update visualization)
            notifyEnvUpdated(this);
        end
        
        % Reset environment to initial state and output initial observation
        function InitialObservation = reset(this)
            % Theta (+- .05 rad)
            T0 = 2 * 0.05 * rand - 0.05;  
            % Thetadot
            Td0 = 0;
            % X 
            X0 = 0;
            % Xdot
            Xd0 = 0;
            
            InitialObservation = [T0;Td0;X0;Xd0];
            this.State = InitialObservation;
            
            % (optional) use notifyEnvUpdated to signal that the 
            % environment has been updated (e.g. to update visualization)
            notifyEnvUpdated(this);
        end
    end
    %% Optional Methods (set methods' attributes accordingly)
    methods               
        % Helper methods to create the environment
        % Discrete force 1 or 2
        function force = getForce(this,action)
            if ~ismember(action,this.ActionInfo.Elements)
                error('Action must be %g for going left and %g for going right.',-this.MaxForce,this.MaxForce);
            end
            force = action;           
        end
        % update the action info based on max force
        function updateActionInfo(this)
            this.ActionInfo.Elements = this.MaxForce*[-1 1];
        end
        
        % Reward function
        function Reward = getReward(this)
            if ~this.IsDone
                Reward = this.RewardForNotFalling;
            else
                Reward = this.PenaltyForFalling;
            end          
        end
        
        % (optional) Visualization method
        function plot(this)
            % Initiate the visualization
            
            % Update the visualization
            envUpdatedCallback(this)
        end
        
        % (optional) Properties validation through set methods
        function set.State(this,state)
            validateattributes(state,{'numeric'},{'finite','real','vector','numel',4},'','State');
            this.State = double(state(:));
            notifyEnvUpdated(this);
        end
        function set.HalfPoleLength(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','HalfPoleLength');
            this.HalfPoleLength = val;
            notifyEnvUpdated(this);
        end
        function set.Gravity(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','Gravity');
            this.Gravity = val;
        end
        function set.CartMass(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','CartMass');
            this.CartMass = val;
        end
        function set.PoleMass(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','PoleMass');
            this.PoleMass = val;
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
        function set.AngleThreshold(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','AngleThreshold');
            this.AngleThreshold = val;
        end
        function set.DisplacementThreshold(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','DisplacementThreshold');
            this.DisplacementThreshold = val;
        end
        function set.RewardForNotFalling(this,val)
            validateattributes(val,{'numeric'},{'real','finite','scalar'},'','RewardForNotFalling');
            this.RewardForNotFalling = val;
        end
        function set.PenaltyForFalling(this,val)
            validateattributes(val,{'numeric'},{'real','finite','scalar'},'','PenaltyForFalling');
            this.PenaltyForFalling = val;
        end
    end
    
    methods (Access = protected)
        % (optional) update visualization everytime the environment is updated 
        % (notifyEnvUpdated is called)
        function envUpdatedCallback(this)
        end
    end
end