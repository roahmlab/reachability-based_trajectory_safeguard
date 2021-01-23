classdef EnvironmentTransform < rl.env.MATLABEnvironment
    % EnvironmentTransform: Transform the observation, action or reward of
    % an environment
    
    % Copyright 2019 The MathWorks, Inc.    
        
    properties (Access = protected)
        % Environment being transformed
        Environment

    end
   
    methods  
        
        function this = EnvironmentTransform(Env)
            % Create a reinforcemnet learning environment for OpenAI Gym
            % OpenAIGymInterface is an OpenAIGymInterface class
            % ObservationInfo is the standard observation description
            % ActionInfo is the standard action description
            
            this = this@rl.env.MATLABEnvironment(Env.ObservationInfo,Env.ActionInfo);
            this.Environment = Env;


        end
        
        function [Observation,Reward,isTerminal,info] = step(this, action)
            error('step must be implemented')

        end  
        
        function InitialState = reset(this)
            error('reset must be implemented')
        end       
          
    end
end