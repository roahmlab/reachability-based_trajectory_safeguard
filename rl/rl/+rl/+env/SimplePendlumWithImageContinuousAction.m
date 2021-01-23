classdef SimplePendlumWithImageContinuousAction < rl.env.AbstractSimplePendlumWithImage
% SIMPLEPENDLUMWITHIMAGECONTINUOUSACTION

% Revised: 11-8-2018
% Copyright 2018 The MathWorks, Inc.

    methods
        function this = SimplePendlumWithImageContinuousAction()
            actionInfo = rlNumericSpec([1 1],...
                'LowerLimit',-2,'UpperLimit',2);
            actionInfo.Name = 'torque';
            this = this@rl.env.AbstractSimplePendlumWithImage(actionInfo);
        end
    end
    methods (Access = protected)
        function setMaxTorque_(this,val)
            this.ActionInfo.UpperLimit =  val;
            this.ActionInfo.LowerLimit = -val;
        end
    end
end