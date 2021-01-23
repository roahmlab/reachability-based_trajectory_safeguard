classdef AbstractSimplePendlumWithImage < rl.env.MATLABEnvironment
% ABSTRACTSIMPLEPENDLUMWITHIMAGE: Creates abstract class for a Simple
% Pendulum model with an image observation. Implements the simple pendulum
% dynmaics:
%
%   (I + m*l^2)*ddtheta = tau - c*dtheta + m*g*l*sin(theta)

% Revised: 11-8-2018
% Copyright 2018 The MathWorks, Inc.

    properties
        Mass = 1
        
        RodLength = 1
        
        RodInertia = 0
        
        Gravity = 9.81
        
        DampingRatio = 0
        
        MaximumTorque = 2
        
        % sample time
        Ts = 0.05
        
        State = zeros(2,1)
        
        % Reward Weights
        Q = diag([1 0.1])
        R = 1e-3
    end
    properties (Transient,Access = private)
        Visualizer = []
    end
    methods (Abstract,Access = protected)
        setMaxTorque_(this,val)
    end
    methods
        function this = AbstractSimplePendlumWithImage(ActionInfo)
            % Define observation info
            ObservationInfo(1) = rlNumericSpec([50 50 1],'LowerLimit',0,'UpperLimit',1);
            ObservationInfo(1).Name = 'pendImage';
            ObservationInfo(2) = rlNumericSpec([1 1]);
            ObservationInfo(2).Name = 'angularRate';
            this = this@rl.env.MATLABEnvironment(ObservationInfo,ActionInfo);
            this.MaximumTorque = 2;
        end
        function set.State(this,state)
            validateattributes(state,{'numeric'},{'finite','real','vector','numel',2},'','State');
            this.State = state(:);
            notifyEnvUpdated(this);
        end
        function set.MaximumTorque(this,val)
            validateattributes(val,{'numeric'},{'real','positive','scalar'},'','MaximumTorque');
            this.MaximumTorque = val;
            setMaxTorque_(this,val);
        end
        function set.Mass(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','Mass');
            this.Mass = val;
        end
        function set.RodLength(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','RodLength');
            this.RodLength = val;
        end
        function set.Gravity(this,val)
            validateattributes(val,{'numeric'},{'finite','real','scalar'},'','Gravity');
            this.Gravity = val;
        end
        function set.DampingRatio(this,val)
            validateattributes(val,{'numeric'},{'finite','real','scalar'},'','DampingRatio');
            this.DampingRatio = val;
        end
        function set.Q(this,val)
            validateattributes(val,{'numeric'},{'real','finite','size',[2 2]},'','Q');
            this.Q = val;
        end
        function set.R(this,val)
            validateattributes(val,{'numeric'},{'real','finite','size',[1 1]},'','R');
            this.R = val;
        end
        function set.Ts(this,val)
            validateattributes(val,{'numeric'},{'finite','real','positive','scalar'},'','Ts');
            this.Ts = val;
        end
        function varargout = plot(this)
            if isempty(this.Visualizer) || ~isvalid(this.Visualizer)
                this.Visualizer = rl.env.viz.SimplePendlumWithImageVisualizer(this);
            else
                bringToFront(this.Visualizer);
            end
            if nargout
                varargout{1} = this.Visualizer;
            end
        end
        function [nextobs,rwd,isdone,loggedSignals] = step(this,action)
            
            loggedSignals = [];
            
            % trapezoidal integration
            ts = this.Ts;
            x1 = this.State;
            [dx1,tau] = dynamics(this,x1,action);
            dx2 = dynamics(this,x1 + ts*dx1,action);
            x = ts/2*(dx1 + dx2) + x1;
            
            % set the state
            this.State = x;
            
            % wrap the output angle
            x(1) = atan2(sin(x(1)),cos(x(1)));
            
            % generate the output image
            img = generateImage(this);
            
            % assign to the observations
            nextobs = {img,x(2)};
            
            % calculate reward
            rwd = - x'*this.Q*x - tau'*this.R*tau;
            
            % terminate at max episodes
            isdone = false;
        end
        function obs = reset(this)
            x = [pi;0];
            this.State = x;
            obs = {generateImage(this),x(2)};
        end
        function img = generateImage(this)
            
            x = this.State;
            
            oinfo = getObservationInfo(this);
            sz = oinfo.Dimension;
            center = ceil(sz(1)/2);
            L = center/2;
            
            s = sin(x(1));
            c = cos(x(1));
            x_ = -round(L*s) + center + (-3:3);
            y_ = -round(L*c) + center + (-3:3);
            
            img = ones(sz);
            img(y_,x_) = 0;
        end
    end
    methods (Access = private)
        function [dx,tau] = dynamics(this,x,tau)
            theta  = x(1);
            dtheta = x(2);
            tau = max(-this.MaximumTorque,min(this.MaximumTorque,tau));
            
            I = this.RodInertia;
            L = this.RodLength;
            m = this.Mass;
            c = this.DampingRatio;
            g = this.Gravity;
            
            I_ = I + m*L^2;
            
            dx(1,1) = dtheta;
            dx(2,1) = 1/I_*(tau - c*dtheta + m*g*L*sin(theta));
        end
    end
end