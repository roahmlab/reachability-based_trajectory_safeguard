classdef CartPoleContinuousAction < rl.env.CartPoleAbstract
    % RLENVCARTPOLECONTINUOUSACTION: Creates cartpole RL environment with 
    % continuous action -Force and Force. Reward for not falling is 1 and
    % penalty for falling is -50 which are different from cartpole
    % environment with discrete actions.
    
    % Copyright 2017-2018 The MathWorks Inc.    
    
    methods        
        function this = CartPoleContinuousAction()
            ActionInfo = rlNumericSpec([1 1]);
            ActionInfo.Name = 'CartPole Action';
            this = this@rl.env.CartPoleAbstract(ActionInfo);
            
            this.RewardForNotFalling = 1;
            this.PenaltyForFalling = -50;
            
            updateActionInfo(this);
        end      
    end
    
    methods (Access = protected)
        % Continuous force between [-this.Force,this.Force]
        function Force = getForce(this,action)
            Force = max(min(action,this.MaxForce),-this.MaxForce);
        end
        % update the action info based on max force
        function updateActionInfo(this)
            this.ActionInfo.LowerLimit = -this.MaxForce;
            this.ActionInfo.UpperLimit =  this.MaxForce;
        end
        % Further from the goal less reward        
        function Reward = getReward(this,x,~)
            distReward = 1 - abs(x)/this.XThreshold;
            if ~this.IsDone
                Reward = 0.5 * this.RewardForNotFalling + 0.5 * distReward;
            else
                Reward = this.PenaltyForFalling;
            end            
        end        
    end    
end