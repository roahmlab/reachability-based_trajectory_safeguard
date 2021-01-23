classdef SimplePendlumWithImageDiscreteAction < rl.env.AbstractSimplePendlumWithImage
% SIMPLEPENDLUMWITHIMAGEDISCRETEACTION

% Revised: 11-8-2018
% Copyright 2018 The MathWorks, Inc.

    methods
        function this = SimplePendlumWithImageDiscreteAction()
            actionInfo = rlFiniteSetSpec([-2 0 2]);
            actionInfo.Name = 'torque';
            this = this@rl.env.AbstractSimplePendlumWithImage(actionInfo);
        end
    end
    methods (Access = protected)
        function setMaxTorque_(this,val)
            this.ActionInfo.Elements = [-1 -0.5 0 0.5 1]*val;
        end
    end
end