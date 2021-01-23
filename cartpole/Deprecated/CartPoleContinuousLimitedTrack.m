% classdef CartPoleContinuousLimitedTrack < rl.env.CartPoleAbstract
classdef CartPoleContinuousLimitedTrack < CartPoleAbstractLimitedTrack
    % RLENVCartPoleContinuousLimitedTrack: Creates cartpole RL environment with 
    % continuous action -Force and Force. Reward for not falling is 1 and
    % penalty for falling is -50 which are different from cartpole
    % environment with discrete actions. 
    % Limited track means we add constraint on x position
    
    % Copyright 2017-2018 The MathWorks Inc.    
    
    methods        
        function this = CartPoleContinuousLimitedTrack()
            ActionInfo = rlNumericSpec([1 1]);
            ActionInfo.Name = 'CartPole Action';
            this = this@CartPoleAbstractLimitedTrack(ActionInfo);
%             this = this@rl.env.CartPoleAbstract(ActionInfo);

%             this.RewardForNotFalling = 1;   % maybe we could change the reward
%             this.PenaltyForFalling = -50;
            
            updateActionInfo(this);
        end      
    end
    
    methods (Access = protected)
        % Continuous force between [-this.Force,this.Force]
        function Force = getForce(this,action)
%             Force = max(min(action,this.MaxForce),-this.MaxForce);
            
            % here we introduce a LLC as the force
            k1 = this.k1;
            k2 = this.k2;
            state = this.State;
%             x = state(1);
            x_dot = state(2);
%             Force = k1*(action-x) - k2*x_dot;
            Force = k1*action - k2*x_dot;
        end
        % update the action info based on max force
        function updateActionInfo(this)
            this.ActionInfo.LowerLimit = -this.MaxForce;
            this.ActionInfo.UpperLimit =  this.MaxForce;
        end
        % Further from the goal less reward        
        function Reward = getReward(this,x,force)
%             distReward = 1 - abs(x)/this.XThreshold;
%             if ~this.IsDone
%                 Reward = 0.5 * this.RewardForNotFalling + 0.5 * distReward;
%             else
%                 Reward = this.PenaltyForFalling;
%             end
            % we use a new reward function
            if this.IsDone
                B = 1;
            else
                B = 0;
            end
            theta = this.State(3);
%             % since we need theat to be the multiple of 2*pi, we cannot
%             % just use theta in the reward function, but the distance of
%             % theta to 2k*pi
%             k = ceil(theta/2/pi);
%             theta = min(theta-2*pi*(k-1),2*pi*k-theta); % this is the distance of theta to the nearest 2*k*pi
            Reward = -0.1*(5*theta^2+x^2+0.05*force^2)-100*B;
        end        
    end    
end