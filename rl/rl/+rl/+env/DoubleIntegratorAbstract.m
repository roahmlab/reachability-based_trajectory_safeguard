classdef DoubleIntegratorAbstract < rl.env.MATLABEnvironment    
    % RLENVDOUBLEINTEGRATORABSTRACT: Creates abstract class for double
    % integrator RL environment.
    
    % Copyright 2017-2018 The MathWorks Inc.
    
    properties
        % gain for the double integrator
        Gain = 1.0
        
        % sample time
        Ts = 0.1
        
        % Max distance, will terminate if absolute position is exceeded
        MaxDistance = 5
        
        % Goal threshold, will terminate if norm of state is less than
        % threshold
        GoalThreshold = 1e-2
        
        % Reward Weights (in continuous time to remove dependence on sample
        % time)
        % reward =  - integral(x'*Q*x + u'*R*u)
        Q  = diag([10 1]) %infinite horizon Q, used to find Qd
        R  = 0.01
    end
    properties (Access = protected)
        MaxForce_ = Inf
    end
    properties (Dependent)
        % Max force
        MaxForce
    end
    properties 
        % system state [x,dx]'
        State = zeros(2,1)
    end
    properties (Transient,Access = private)
        Visualizer = []
    end
    properties (Access = private)
        % reward = - sum_i(x_i'*Qd*x_i + u_i'*Rd*u_i)
        % Note, these weights are derived from Q and R and will be updated
        % after Q and R are set
        Qd (2,2) double = eye(2)
        Rd (1,1) double = 0.01
        Nd (2,1) double = [0 0]'
        
        % discrete system matrices
        Ad (2,2) double
        Bd (2,1) double
        Cd (2,2) double        
    end
    methods (Abstract, Access = protected)
        force = getForce(this,force);
    end
    
    methods (Access = protected)
        function setMaxForce_(this,val)
            % define how the setting of max force will behave, which is
            % different among continuous/discrete implementations of the
            % environment
            this.MaxForce_ = val;
            this.ActionInfo.Values = [-val,val];
        end
        function updatePerformanceWeights(this)            
            % get the continuous linearized system
            a = [0 1;0 0];
            b = [0;this.Gain];
            
            % Determine discrete equivalent of continuous cost function
            % along with Ad,Bd matrices of discretized system. This is
            % equivalent to the discretization in lqrd
            q = this.Q;
            r = this.R;
            nn = [0 0]';
            Nx = 2; Nu = 1;
            n = Nx+Nu;
            Za = zeros(Nx); Zb = zeros(Nx,Nu); Zu = zeros(Nu);
            M = [ -a' Zb   q  nn
                -b' Zu  nn'  r
                Za  Zb   a   b
                Zb' Zu  Zb' Zu];
            phi = expm(M*this.Ts);
            phi12 = phi(1:n,n+1:2*n);
            phi22 = phi(n+1:2*n,n+1:2*n);
            QQ = phi22'*phi12;
            QQ = (QQ+QQ')/2;        % Make sure QQ is symmetric
            qd = QQ(1:Nx,1:Nx);
            rd = QQ(Nx+1:n,Nx+1:n);
            nd = QQ(1:Nx,Nx+1:n);
            ad = phi22(1:Nx,1:Nx);
            bd = phi22(1:Nx,Nx+1:n);
            
            this.Rd = rd;
            this.Qd = qd;
            
            this.Ad = ad;
            this.Bd = bd;
            this.Nd = nd;
            this.Cd = eye(2);
        end        
    end
    methods
        function this = DoubleIntegratorAbstract(ActionInfo)
            % Define observation info
            ObservationInfo = rlNumericSpec([2 1]);
            ObservationInfo.Name = 'states';
            ObservationInfo.Description = 'x, dx';
            this = this@rl.env.MATLABEnvironment(ObservationInfo,ActionInfo);
            updatePerformanceWeights(this);
        end
        function set.State(this,state)
            validateattributes(state,{'numeric'},{'finite','real','vector','numel',2},'','State');
            this.State = state(:);
            notifyEnvUpdated(this);
        end
        function set.GoalThreshold(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','GoalThreshold');
            this.GoalThreshold = val;
        end
        function set.MaxDistance(this,d)
            validateattributes(d,{'numeric'},{'finite','real','positive','scalar'},'','MaxDistance');
            this.MaxDistance = d;
            notifyEnvUpdated(this);
        end
        function varargout = plot(this)
            if isempty(this.Visualizer) || ~isvalid(this.Visualizer)
                this.Visualizer = rl.env.viz.DoubleIntegratorVisualizer(this);
            else
                bringToFront(this.Visualizer);
            end
            if nargout
                varargout{1} = this.Visualizer;
            end
        end
        function set.MaxForce(this,val)
            validateattributes(val,{'numeric'},{'real','positive','scalar'},'','MaxForce');
            setMaxForce_(this,val);
        end
        function val = get.MaxForce(this)
            val = this.MaxForce_;
        end
        function set.Q(this,val)
            validateattributes(val,{'numeric'},{'real','finite','size',[2 2]},'','Q');
            this.Q = val;
            updatePerformanceWeights(this);
        end
        function set.R(this,val)
            validateattributes(val,{'numeric'},{'real','finite','size',[1 1]},'','R');
            this.R = val;
            updatePerformanceWeights(this);
        end
        function set.Gain(this,val)
            validateattributes(val,{'numeric'},{'finite','real','scalar'},'','Gain');
            this.Gain = val;
            updatePerformanceWeights(this);
        end
        function set.Ts(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','Ts');
            this.Ts = val;
            updatePerformanceWeights(this);
        end
        function [nextobs,rwd,isdone,loggedSignals] = step(this,action)
                       
            loggedSignals = [];
            
            % saturate the action
            action = getForce(this,action);
            
            % get next state
            x = this.State;
            xk1 = this.Ad*x + this.Bd*action;
            nextobs = this.Cd*xk1;
            this.State = xk1;
            
            % calculate reward
            rwd = - x'*this.Qd*x - action'*this.Rd*action - 2*x'*this.Nd*action;
            
            % The episode will terminate under the following conditions:
            % 1. the mass moves more than X units away from the origin
            % 2. the norm of the state is less than some threshold
            %
            % The second point is critical for training as it prevents the
            % replay buffer being saturated with 0s for training
            isdone = abs(x(1)) > this.MaxDistance || norm(x) < this.GoalThreshold;
        end
    end
end