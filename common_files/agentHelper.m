classdef agentHelper < handle
    %% properties
    properties
        A
        % FRS object
        zono_full
        
        %time info
        t_move
        t_failsafe_move
        
        %saved TUZ for failsafe action when none can be found
        T
        U
        Z
        
        % flags to pass into children AH
        flags
        
        stuck_count = 0
        verbose = 0 ;
        plot_data
        name = 'agentHelper' ;

        time_vec = [];
        timeout_err_counter
    end
    %% methods
    methods
        function AH = agentHelper(A,FRS_path,varargin)
            % AH = agentHelper(A,FRS_path,varargin) Construct by passing in
            % a agent and FRS file, file handled by children class 
            if nargin > 0
                AH = parse_args(AH,varargin{:}) ;
                AH.A = A;
                AH.zono_full = load(FRS_path);%'zono_full_7.13_1spd.mat');
            end
        end
        %% the only useful function in the parent class
        function [action_replaced, replace_distance, stop_sim, k] = advanced_move(AH,action,world_info)
            % [action_replaced, replace_distance, stop_sim, k] = advanced_move(AH,action,world_info)
            % pass in a action from the main file with action range of [-1,1]
            % scale properly and move the agent accordingly, checks for
            % action safety if safety layer is on
            start_plan_tic = tic;
            
            k = AH.convert_action_to_parameter(action,AH.flags.discrete_flag);
            stop_sim = 0;
            if strcmp(AH.flags.safety_layer, 'Z') || strcmp(AH.flags.safety_layer, 'RTS')%zonotope safety mode, do learning and then use the agent
                [k, replace_distance, replaced_flag]= AH.adjust(k,world_info);
                if replaced_flag == 0
                    no_replace_action = 0;
                    action_replaced = 0;
                elseif replaced_flag == 1
                    no_replace_action = 0;
                    action_replaced = 1;
                elseif replaced_flag == 2
                    no_replace_action = 1;
                    action_replaced = 1;
                end
            elseif strcmp(AH.flags.safety_layer, 'N') || strcmp(AH.flags.safety_layer, 'NoSafety')%no safety but still has learning
                replace_distance = 0;
                action_replaced = 0;
                no_replace_action = 0;
            elseif strcmp(AH.flags.safety_layer, 'R') || strcmp(AH.flags.safety_layer, 'RTD')%RTD mode, use any agent, does not use agent output, just run planning
                %uses a high level planner(HLP), to supply waypoint and
                %then optimize in trajecotry space, doesn't need learning.
                %can supply a naive HLP that always output the goal position for a direct comparison with RL. 
                [k, no_replace_action] = AH.gen_param(world_info);
%                 [k_replaced, replace_distance, ~]= AH.adjust(k,world_info);
                replace_distance = 0;
                action_replaced = 0;
            end

            AH.time_vec =[AH.time_vec toc(start_plan_tic)];

            if ~no_replace_action 
                AH.stuck_count = 0;
                [AH.T, AH.U, AH.Z] = AH.gen_ref(k);
                if AH.plot_flag % this is for highway only, may have bug for drone
                    h1 = figure(1); 
                    AH.plot();
                    AH.A.plot();
                    state=AH.A.state(:,end);
%                     text(state(1)+20,2,[num2str(state(4)) 'm/s'],'Color','red','FontSize',15)
                    figure(9);
                    set(gca,'FontSize',15)
                    axl=subplot(1,2,1);
                    axis equal
                    set(gca,'FontSize',15)
                    copyobj(h1.Children.Children,axl)
                    axl.XLim = h1.Children.XLim;axl.YLim = h1.Children.YLim;
                    if axl.XLim(1)+15 > axl.XLim(2)-30
                        axl.XLim = [axl.XLim(1) axl.XLim(2)];
                    else
                        axl.XLim = [axl.XLim(1)+15 axl.XLim(2)-30];
                    end
                        
                    xlabel('x [m]');
                    ylabel('y [m]');
                    
                end
                AH.A.move(AH.t_move,AH.T, AH.U, AH.Z);
            else
                fprintf('no replacement action found\n');
                AH.stuck_count = AH.stuck_count + 1;
                if AH.stuck_count == 1
                    idx = find((AH.T-AH.t_move)>0);
                    idx = idx(1)-1;
                    AH.A.move(AH.t_failsafe_move,AH.T(idx:end)-AH.T(idx), AH.U(:,idx:end), AH.Z(:,idx:end));
                elseif AH.stuck_count < 4
                    AH.T = [0 AH.t_failsafe_move];
                    AH.U = repmat(AH.U(:,end),1,2);
                    AH.Z = repmat(AH.Z(:,end),1,2);
                    AH.A.move(AH.t_failsafe_move,AH.T, AH.U, AH.Z);
                else
                    stop_sim = 1;
                    fprintf('Terminating\n');
                end
            end
        end
        %% Helper and unimplemented functions
        function agent_info = get_agent_info(AH)
            agent_info = AH.A.get_agent_info();
        end
        function reset(AH,flags,world_start)
            if nargin < 3
                AH.A.reset();
            else
                AH.A.reset(world_start);
            end
            AH.flags = flags;
            error('agent helper reset not implemented')
        end % following functions should be implemented in the child class
        function K = convert_action_to_parameter(AH,action)
            error('convert_action_to_parameter not implemented')
        end
        function [K, replace_distance, replaced_flag] = adjust(AH,K,flags)
            error('adjust not implemented')
        end
        function [T, U, Z] = gen_ref(AH,K)
            error('gen_ref not implemented')
        end
        function [k, no_action_found] = gen_param(AH,world_info)
            error('gen_param not implemented')
        end
        function  plot(AH)
            error('plot not implemented')
        end
        function plot_adjust(AH)
            error('plot_adjust not implemented')
        end
        function plot_A(AH)
            AH.A.plot()
        end

        function vdisp(AH,s,l)
            % Display a string s if the message's verbose level l is greater
            % than or equal to the planner's verbose level.
            if nargin < 3
                l = 1 ;
            end
            if AH.verbose >= l
                if ischar(s)
                    disp(['    AH: ',s])
                else
                    disp('    AH: String not provided!')
                end
            end
        end
    end
end
