classdef rlsimulator < handle
    %% properties
    properties
        safety_layer %N: none, Z: zonotope safety
        discrete_flag = 0;
        replace_action = 0; % 0 for original rl lib, 1 for replace with new
        %action, 2 for with punishing reward b/w replaced and original
        plot_sim_flag = 1; %plot simulation figure 1
        plot_AH_flag;
        plot_adjust_flag = 0;% plot replacement process figure 2
        eval = 0; % used to set random seed, using the same seed generate same random environment
        AH
        W
        
        fig_num = 1
        epscur = 1;
        time_vec=[];
        iter_num = 0;
    end
    %% methods
    methods
        %% constructor
        function S = rlsimulator(AH,W,varargin)
            S = parse_args(S,varargin{:}) ;
            S.AH = AH;
            S.W  = W;
            if isempty(S.safety_layer)
                error('Please specify to use safety layer: Z or not: N')
            end
            if S.replace_action == 0
                rmpath(genpath('./Safe_RL_Highway_SS/rl'));
            elseif S.replace_action == 1
            elseif S.replace_action == 2
            else
                error('Not valid replace action flag')
            end
            figure('units','normalized','outerposition',[0 0 1 1])
        end
        
        %% step
        function [Observation,Reward,IsDone,LoggedSignals,varargout] = step(S,action)
            %step in time and check for done flags
            tic
            agent_info = S.AH.get_agent_info();
            world_info = S.W.get_world_info(agent_info);%where obstacles are
            % move
           
            [action_replaced, replace_distance, stuck, k] = S.AH.advanced_move(action,world_info);
            agent_info = S.AH.get_agent_info() ;
            agent_info.replace_distance = replace_distance;
            Observation = S.W.get_ob(agent_info);
            Reward = S.W.getRew(agent_info,Observation);
            collision = S.W.collision_check(agent_info);
            goal_check = S.W.goal_check(agent_info);
            IsDone = S.determine_isDone_flag(collision,action_replaced,stuck,goal_check);
            if IsDone == 5
                Reward  = Reward +100;
            end
            if S.eval &&( IsDone == 1 || IsDone == 3 ||IsDone == 4 || IsDone == 5) 
               Filename = sprintf('sim_summary_%s-%s.mat', [num2str(IsDone) strcat(num2str(S.epscur),"_",datestr(now,'HH-MM-SS.FFF'))]);
               % reference data for plot
               size(S.time_vec)
               ref_Z = S.AH.ref_Z;
               proposed_ref_Z = S.AH.proposed_ref_Z;
               T= S.AH.T;
               t_real_start_arr = S.AH.t_real_start;
               t_proposed_start_arr = S.AH.t_proposed_start;
               t_move =  S.AH.t_move;
               plotting_param = S.AH.FRS_plotting_param;
               save(Filename,'agent_info','world_info','ref_Z','proposed_ref_Z','plotting_param','T','t_move','t_real_start_arr','t_proposed_start_arr')
            end
             S.plot();
            LoggedSignals = struct;
%             LoggedSignals.secondexp = single(zeros(size(action,1),1)); %this will be ignored by original rl library
            %but will be pickup up by edited rl library to train the agent
            %differently. Here just dummy to make it not error.

            % send output action if nargout > 4
            if nargout > 4
                varargout = {action_replaced,k} ;
            end
            t_step = toc;
            S.time_vec = [S.time_vec t_step];
            
        end
        
        %% reset
        function [iniOb, LoggedSignals] = reset(S)
            % reset if eval different
            if S.eval
            S.epscur
            end
            LoggedSignals = struct;
            flags = struct;
            flags.discrete_flag = S.discrete_flag;
            flags.replace_action = S.replace_action;
            flags.safety_layer = S.safety_layer;
            if S.eval 
               S.W.setup(S.epscur);
               S.AH.reset(flags,S.epscur);
            else
               S.W.setup();
               S.AH.reset(flags);% inconsistent between eval and not, fix eval for drone on this
            end
            
            iniOb = S.W.get_ob(S.AH.get_agent_info());
            S.plot();
            S.epscur = S.epscur + 1;
            S.time_vec = [];
        end
        
        %% helper functions
        function plot(S)
            % plot
            if S.plot_sim_flag 
                fh = figure(S.fig_num) ;
                cla ; hold on ; axis equal ;
                S.W.plot();
                
                xlabel('x [m]');
                ylabel('y [m]');
                set(gca,'FontSize',15)

                if  S.plot_AH_flag
                    S.AH.plot()
                    if S.plot_adjust_flag 
                    S.AH.plot_adjust()
                    end
                end
                S.AH.plot_A();
            end
            
            
            
        end
       
        function [IsDone] = determine_isDone_flag(S,collision,action_replaced,stuck,goal_check)
            if collision && action_replaced
                IsDone = 3;
            elseif collision
                IsDone = 1;
            elseif goal_check
                IsDone = 5;
            elseif stuck
                IsDone = 4;
            elseif action_replaced
                IsDone = 2;
            else
                IsDone = 0;
            end  
        end
    end
end
