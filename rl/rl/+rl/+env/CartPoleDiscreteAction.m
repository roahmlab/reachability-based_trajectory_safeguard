classdef CartPoleDiscreteAction < rl.env.CartPoleAbstract
    % RLENVCARTPOLEDISCRETEACTION: Creates cartpole RL environment with 
    % discrete action (1,left or 2, right). Reward for not falling is 1 and
    % penalty for falling is -5 which are different from cartpole
    % environment with continuous actions.
    
    % Copyright 2017-2018 The MathWorks Inc.    
        
    methods
        function this = CartPoleDiscreteAction()
            ActionInfo = rlFiniteSetSpec([-1 1]);
            ActionInfo.Name = 'CartPole Action';
            this = this@rl.env.CartPoleAbstract(ActionInfo);
            
            updateActionInfo(this);
        end        
    end
    
    methods (Access = protected)
        % Discrete force 1 or 2
        function force = getForce(this,action)
            if ~ismember(action,this.ActionInfo.Elements)
                error(message('rl:env:CartPoleDiscreteInvalidAction',sprintf('%g',-this.MaxForce),sprintf('%g',this.MaxForce)));
            end
            force = action;           
        end
        % update the action info based on max force
        function updateActionInfo(this)
            this.ActionInfo.Elements = this.MaxForce*[-1 1];
        end
        function Reward = getReward(this,~,~)
            if ~this.IsDone
                Reward = this.RewardForNotFalling;
            else
                Reward = this.PenaltyForFalling;
            end          
        end        
    end
end