classdef RewardTransform < rl.env.EnvironmentTransform
    % RewardTransform: Transform the reward of an environment
    
    % Copyright 2019 The MathWorks Inc.

    methods
        function [Observation,Reward,isTerminal,info] = step(this, action)
            % Revisit: need to consider upstream call from nargin
            [Observation,Reward,isTerminal,info] = step(this.Environment, action);
            Reward = transform(this,Reward);
        end
        
        function InitialState = reset(this)
            InitialState = reset(this.Environment);
        end
    end
    
    methods (Access = protected)
        function this = RewardTransform(env)
            this = this@rl.env.EnvironmentTransform(env);        
        end
    end
    
    methods (Abstract)
        Reward = transform(this,Reward)       
    end

end