classdef cartpole_sim_agent > RTD_agent_2D
    properties
        integrator_type = 'ode45' ; % choose 'ode45' or 'ode4' or 'ode113'
        integrator_time_discretization = 0.01 ; % for ode4
        desired_initial_condition=[-490; 6; 0; 22; 0];
        
    end
    methods
        function A = cartpole_sim_agent(varargin)
            name = 'highway_cruiser' ;
            default_footprint = [4.8 2];
            n_states = 3 ;
            n_inputs = 1 ;
            A@RTD_agent_2D('name',name,...
                'footprint',default_footprint,...
                'n_states',n_states,'n_inputs',n_inputs,...
                'stopping_time',stopping_time,'sensor_radius',sensor_radius,varargin{:}) ;