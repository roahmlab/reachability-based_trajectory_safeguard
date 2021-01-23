classdef ActionTransform < rl.env.EnvironmentTransform
    % ActionTransform: Transform the action of an environment
    
    % Copyright 2019 The MathWorks Inc.

    methods
        function [Observation,Reward,isTerminal,info] = step(this, action)
            % Revisit: need to consider upstream call from nargin
            [Observation,Reward,isTerminal,info] = step(this.Environment, transform(action));
        end
        
        function InitialState = reset(this)
            InitialState = reset(this.Environment);
        end
    end
    
    methods (Access = protected)
        function this = ActionTransform(env)
            this = this@rl.env.EnvironmentTransform(env);        
        end
    end
    
    methods (Abstract)
        ActionTransform = transform(this,ActionTransform)       
    end

end