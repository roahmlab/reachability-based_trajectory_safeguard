classdef DoubleIntegratorDiscreteAction < rl.env.DoubleIntegratorAbstract
    % RLENVDOUBLEINTEGRATORDISCRETEACTION: Creates double integrator RL
    % environment with discrete action.

    % Copyright 2017-2018 The MathWorks Inc.
    
    methods (Access = protected)
        function force = getForce(~,force)
            % saturate the action
        end
        
        function setMaxForce_(this,val)
            validateattributes(val,{'numeric'},{'scalar','finite'});
            this.MaxForce_ = val;
            this.ActionInfo.Elements = [-val,0,val];
        end
    end
    
    methods
        function this = DoubleIntegratorDiscreteAction()
            ActionInfo = rlFiniteSetSpec([-2 0 2]);
            ActionInfo.Name = 'force';
            this = this@rl.env.DoubleIntegratorAbstract(ActionInfo);
            this.MaxForce = 2;
        end        

        function obs = reset(this)
            % reset the double integrator to +/- 0.4*MaxDistance
            r = rand(1,1) - 0.5;
            if r == 0, r = 1; end
            % X
            X0 = 0.4*this.MaxDistance*sign(r);
            % Xdot
            Xd0 = 0;
            this.State = [X0;Xd0];
            obs = this.State;
        end
    end
end