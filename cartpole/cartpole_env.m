classdef cartpole_env < handle
    properties
        Gravity = 9.8;
        MassCart = 1.0;
        MassPole = 0.1;
        Length = 0.5;
        MaxForce = 1;
        Ts = 0.1;
        Tf = 2;
        ThetaThresholdRadians = 12 * pi/180; % 12 * pi/180
        XThreshold = 3.5;    %2.4;
        max_force = -inf;
        max_x_acc = -inf;
        plottor_value = 0;
        
        plottor;
        plottor_flag = 1;
       
        
        % Define the error bound
        error_bound_min = [];
        error_bound_max = [];
        
        % define the sliced resolution and range
        % set the number of points we want to sample
        N_x = 5;
        N_x_dot = 10;
        N_theta = 10;
        N_theta_dot =10;
        % set the range of sampled space
        x_start = -2; x_end = 2;
        x_dot_start = -5; x_dot_end = 5;
        theta_start = 0; theta_end = 2*pi;
        theta_dot_start = -5; theta_dot_end = 5;
        finish_count = 0;
        
        evalflag = 0;
        eval_iter = 0;

        AH
        
    end
    methods
        function this = cartpole_env(AH)
            this.AH = AH;
%             cartpole_agentHelper(Env,A,'cartpole_frs_9_14.mat','safety_layer','N');
        end
        
        
        function initialobservation = reset(this,ic)
            this.eval_iter  %%%%%%%%% print
            % Randomize the initial pendulum angle between (+- .05 rad)
            
            % Theta (+- .05 rad)
            %             T0 = 2*0.05*rand - 0.05;
            % in order to swing up, we put the initial condition as pi
          
            
            if this.evalflag
                this.eval_iter = this.eval_iter +1;
                rng(this.eval_iter);
                this.AH.A.reset([0;0;-pi+2*pi*rand(); -1+2*rand()]);
            else
                this.AH.A.reset();
            end
            agent_info = this.AH.A.get_agent_info();
            state = agent_info.state(:,end);
            X0 = state(1);
            Xd0 = state(2);
            T0 = state(3);
            Td0 = state(4); 
            %             initialState= [T0;Td0;X0;Xd0];
            %             initialState= [X0;T0;Xd0;Td0];
            initialobservation = [sin(T0);cos(T0);Td0;X0;Xd0];

            
            
        end
        
        function [Observation,Reward,IsDone,LoggedSignals] = step(this,action)
            
            LoggedSignals = [];
            if this.plottor_flag == 1
                if this.plottor_value == 0
                    this.plottor = rl.env.viz.CartPoleVisualizer(this);
                    this.plottor_value = this.plottor_value +1;
                end
            end
            % Get action
%             this.AH.kv = this.v;
%             this.AH.ka = this.a;
            
            world_info = struct;
            world_info.obstacles = [-5 -3 NaN 3 5];
%             world_info.ka = this.a;
            [action_replaced, ~, stop_sim, k]= this.AH.advanced_move(action,world_info);
            agent_info = this.AH.A.get_agent_info();
            state = agent_info.state(:,end);
%             time = agent_info.time;
            
%             this.a = (state(2,end)- state(2,end-1))/(time(end)-time(end-1));

%             this.State = state(:,end);

            Observation(1) = sin(state(3));
            Observation(2) = cos(state(3));
            Observation(3) = state(4);
            Observation(4) = state(1);
            Observation(5) = state(2);
            
            Observation= Observation';
            
            if (abs(Observation(2) - 1) < 0.1) && (abs(Observation(1) - 0) < 0.1)
                this.finish_count = this.finish_count +1;
            else
                this.finish_count = 0;
            end
           
            
            %state order :X0;Xd0;T0;Td0
            %Observation order: sin(T0);cos(T0);Td0;X0;Xd0
            
            x = Observation(4);
            if stop_sim
                IsDone = 4;
                Filename = sprintf('sim_summary_%s-%s.mat', [num2str(IsDone) strcat(num2str(this.eval_iter),"_",datestr(now,'HH-MM-SS.FFF'))]);
                save(Filename);
            elseif action_replaced && (abs(x) > this.XThreshold)
                IsDone = 3;
                Filename = sprintf('sim_summary_%s-%s.mat', [num2str(IsDone) strcat(num2str(this.eval_iter),"_",datestr(now,'HH-MM-SS.FFF'))]);
                save(Filename);
            elseif ~action_replaced && (abs(x) > this.XThreshold)
                IsDone = 1;
                Filename = sprintf('sim_summary_%s-%s.mat', [num2str(IsDone) strcat(num2str(this.eval_iter),"_",datestr(now,'HH-MM-SS.FFF'))]);
                save(Filename);
            elseif action_replaced && (abs(x) < this.XThreshold)
                IsDone = 2;
                LoggedSignals.secondexp = k;
            elseif this.finish_count == 100
                IsDone = 5;
%                 IsDone = 5;
            else
                IsDone = 0;
            end
            
            if this.plottor_flag == 1
                this.plottor.plot();
            end
            
            % Get reward
            Reward = this.getReward(state,IsDone);
        end
        
        function Force = getForce(this,action)
            Force = max(min(action,this.MaxForce),-this.MaxForce);
        end
        
        function Reward = getReward(this,state,IsDone)
            
            if IsDone == 1
                B = 1
            else
                B = 0;
            end
            theta = state(3);
            
            Reward = (cos(theta)+1)/2 - 0.1*sign(state(1))*sign(state(2)) -30*B-0.05*abs(state(4));% - 0.1* abs(state(1));%-0.01* abs(this.State(4)); %+0.05*force^2
        end
        
    end
end
