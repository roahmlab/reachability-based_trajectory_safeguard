classdef ObservationTransform < rl.env.EnvironmentTransform
    % ObservationTransform: Transform the observation of an environment
    
    % Copyright 2019 The MathWorks Inc.

    
    methods
        function [Observation,Reward,isTerminal,info] = step(this, action)
            % Revisit: need to consider upstream call from nargin
            [Observation,Reward,isTerminal,info] = step(this.Environment, action);
            Observation = transform(this,Observation);
        end
        
        function InitialState = reset(this)
            InitialState = transform(reset(this.Environment));
        end
    end
    
    methods (Access = protected)
        function this = ObservationTransform(env)
            this = this@rl.env.EnvironmentTransform(env);        
        end
    end
    
    methods (Abstract)
        obs = transform(this,obs)       
    end

end