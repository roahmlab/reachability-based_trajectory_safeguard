classdef DoubleIntegratorContinuousAction < rl.env.DoubleIntegratorAbstract
    % RLENVDOUBLEINTEGRATORCONTINUOUSACTION: Creates double integrator RL
    % environment with continuous action.
    
    % Copyright 2017-2018 The MathWorks Inc.
    
    
    methods (Access = protected)
        function force = getForce(this,force)
            % saturate the action
            umax = this.MaxForce;
            force = min(max(force,-umax),umax);
        end
        
        function setMaxForce_(this,val)
            validateattributes(val,{'numeric'},{'scalar'});
            this.MaxForce_ = val;
            this.ActionInfo.UpperLimit =  val;
            this.ActionInfo.LowerLimit = -val;
        end
    end
    
    methods
        function this = DoubleIntegratorContinuousAction()
            % ENV = DOUBLEINTEGRATORCONTINUOUSACTION()
            ActionInfo = rlNumericSpec([1 1]);
            ActionInfo.Name = 'force';
            this = this@rl.env.DoubleIntegratorAbstract(ActionInfo);
        end
        function obs = reset(this)
            % reset the double integrator to +/- 0.8*MaxDistance
            r = rand(1,1) - 0.5;
            if r == 0, r = 1; end
            % X
            X0 = 0.8*this.MaxDistance*sign(r);
            % Xdot
            Xd0 = 0;
            this.State = [X0;Xd0];
            obs = this.State;
        end
    end
end