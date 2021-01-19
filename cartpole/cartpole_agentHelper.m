classdef cartpole_agentHelper < agentHelper
    %% properties
    properties
        
        % property to save the current planned trajectory
        current_trajectory
        current_parameter
        
        % how to get initial condition for previous trajectory
        use_agent_for_initial_condition_flag = true ;
        
        % speed and acceleration info
        v_max = 5 ;
        
        %         k_adjusted = [0;0;0];
        
        plot_cart_flag = 0;
        plot_flag = 0;
        eps = 0.1;
        lowerbound_array = [];
        higherbound_array = [];
        HLP = [];
        finish_count =0;
        fminconopt =[];
        optmal_opt = "sample";
        N_sample = 15;
        N_stop = 0;
        N_stop_threshold = 3;
    end
    
    %% methods
    methods
        %% agentHelper required methods start here
        %% constructor
        function AH = cartpole_agentHelper(A,FRS_path,varargin)
            %              AH@agentHelper(A,FRS_path,varargin{:});
            
            if strcmp(varargin{1},'safety_layer')
                safety_layer = varargin{2};
            end
            AH.flags.discrete_flag=0;
            AH.flags.safety_layer= safety_layer;
            AH.A = A;
            AH.zono_full = load(FRS_path);
            AH.t_move = 0.1;
            AH.t_failsafe_move =0.199;
            o = optimoptions('fmincon') ;
            o.OptimalityTolerance = 1e-6 ;
            o.MaxIterations = 100000 ;
            o.MaxFunctionEvaluations = 10000 ;
            o.SpecifyConstraintGradient = true ;
            o.SpecifyObjectiveGradient = isa(AH.HLP,'cartpole_HLP');
            o.Display = 'off';
            o.CheckGradients = false;
            %             end
            AH.fminconopt = o;
            % Offline
        end
        
        %% reset
        function reset(AH,flags,world_start)
            AH.kv = 0;
            %             AH.ka = 0;
            %             % quadrotor_agentHelper.reset(flags,world_start)
            %             %
            %
            %             %
            %             AH.k_adjusted = [0;0;0];
            %             AH.vdisp('Resetting agent',4)
            %             if nargin < 3
            %                 AH.A.reset() ;
            %             else
            %                 AH.A.reset(world_start) ;
            %             end
            %
            %             AH.vdisp('Setting flags',4)
            %             AH.flags = flags ;
        end
        
        %% adjust
        function [K, replace_distance, replaced_flag] = adjust(AH,k_user,world_info)
            %                 k = k_user; %???
            [cur_ka, kv]=AH.get_ka_kv();
            
            agent_state = AH.A.state(:,end);
            O_all = world_info.obstacles;
            %             x = agent_state(1);
            x_dot = agent_state(2);
            theta = agent_state(3);
            theta_dot = agent_state(4);
            
            %use indicies to select right PRS
            [~,ka_idx] =min (abs(AH.zono_full.ka - cur_ka));
            [~,vd_idx]  = min(abs(AH.zono_full.kpeak - (k_user)));%aka k peak idx
            [~,v_idx] = min(abs(AH.zono_full.vi - x_dot));
            %only vd(k_peak) is relative to current velocity
            
            %plus those following idices to select right ERS
            theta_r = rem(theta,2*pi);%find the remainder
            if theta_r < -pi
                theta_r = theta_r + 2*pi;
            elseif theta_r > pi
                theta_r = theta_r - 2*pi;
            end
            [~,theta_idx] = min(abs([AH.zono_full.theta_a 3.1415] - theta_r));
            if theta_idx == 5
                theta_idx = 1;
            end
            %theta has circular property
            
            [~,theta_dot_idx] = min(abs(AH.zono_full.theta_dot_a - theta_dot));
            
            agent_state_extra = [agent_state; cur_ka];
            [K,replaced_flag] = AH.find_replace_action_global(O_all,k_user,agent_state_extra,v_idx,ka_idx,vd_idx,theta_idx,theta_dot_idx);
            %
            if isempty(K)
                replaced_flag = 2;
            end
            if replaced_flag == 1
                replace_distance = abs(k_user - K);
                %                 fprintf('action %.1f replaced to %.1f',k_user,K)
            else
                replace_distance = 0;% this dist is not gonna be used if 0
            end
        end
        
        function [K, replaced] = find_replace_action_global(AH,O,k_user,agent_state,v_idx,ka_idx,vd_idx,theta_idx,theta_dot_idx)
            %for each combination of v_des_index and y_des_idx, try finding
            %the replacment action. start with the current index, if cannot
            %find coninue onto near ones.
            
            %replace idx in 1d is a line, so just try left, try right....
            for replace_idx = 0:3% 3 v desired to look at, center left and right
                if replace_idx == 0
                    kpkidx = vd_idx; % user selected one check first
                elseif replace_idx == 1
                    kpkidx = v_idx;
                elseif replace_idx == 2
                    kpkidx = v_idx - 1;
                elseif replace_idx == 3
                    kpkidx = v_idx + 1;
                end
                
                
                if kpkidx < 1 ||  kpkidx > length(AH.zono_full.kpeak)
                    %outside of selectable zonotopes, skip
                    continue;
                end
                %                 fprintf('proposing v_peak = %f\n',AH.zono_full.kpeak(kpkidx))
                
                [K, inner_replaced_flag] = AH.find_replace_action(O,k_user,agent_state,v_idx,ka_idx,kpkidx,theta_idx,theta_dot_idx);% find a action that works
                if ~isempty(K)
                    if AH.plot_cart_flag
                        AH.plot_selected_zono(O,K,agent_state,v_idx,ka_idx,kpkidx,theta_idx,theta_dot_idx);
                    end
                    [AH.lowerbound_array, AH.higherbound_array] = AH.traj_selected_zono(O,K,agent_state,v_idx,ka_idx,kpkidx,theta_idx,theta_dot_idx);
                    break
                end
                %                 fprintf('safe action not found\n')
            end
            if replace_idx == 0 && inner_replaced_flag == 0
                replaced = 0;
            else
                replaced = 1;
            end
        end
        function plot_selected_zono(AH,O,K,agent_state,v_idx,ka_idx,vd_idx,theta_idx,theta_dot_idx)
            vd_idx_relative = vd_idx-v_idx+2;%convert from abs index to relative idx.
            zono_peak     = AH.zono_full.res{v_idx,ka_idx,vd_idx_relative,1};
            zono_failsafe = AH.zono_full.res{v_idx,ka_idx,vd_idx_relative,2};
            error_data_cur = AH.zono_full.error_table{v_idx,ka_idx,vd_idx_relative,theta_idx,theta_dot_idx};
            zono_one_all = [zono_peak; zono_failsafe];
            n = length(zono_one_all);
            
            figure(1);clf; hold on;
            pos_cart = agent_state(1);
            O_pts = O - pos_cart;
            %             yline(O_pts(1));yline(O_pts(2));yline(O_pts(4));yline(O_pts(5));
            sh_frs = scatter([],[],30);
            sh_prs = scatter([],[],10);
            for t_idx = 1:n
                
                zono_one = zono_one_all{t_idx}{1};
                
                if isempty(zono_one)
                    error("no such zonotope"+num2str(v_ini_idx)+num2str(y_des_idx)+num2str(h_ini_idx)+num2str(del_idx)+num2str(vd_idx));
                end
                %                 y_des = K(2);
                %                 vx_des = K(1);
                
                %                  P.vx_des,P.y_des - agent_state(2),agent_state(3),agent_state(4),agent_state(5)
                kv_soft = bound_values(agent_state(2),AH.zono_full.vi(v_idx)-AH.zono_full.kvg,AH.zono_full.vi(v_idx)+AH.zono_full.kvg);
                ka_soft = bound_values(agent_state(5),AH.zono_full.ka(ka_idx)-AH.zono_full.kag,AH.zono_full.ka(ka_idx)+AH.zono_full.kag);
                zono_one = zonotope_slice(zono_one, [3;4], [kv_soft;ka_soft]);
                %                zono_one = zonotope_slice(zono_one, [4;5], [vi_soft;agent_state(3)]);
                centererr = mean(error_data_cur(t_idx,:));
                gen_err = zeros(5,5); gen_err(1,1) = (error_data_cur(t_idx,2)- error_data_cur(t_idx,1))/2;
                err_zono  = zonotope([[centererr;0;0;0;0], gen_err]);
                
                
                full_zono = zono_one + err_zono;
                frs_zono = zonotope_slice(full_zono, 2, K);
                prs_zono = zonotope_slice(zono_one, 2, K);
                %Slice, Do it now!
                zono_cur_vert = deleteAligned(project(frs_zono, 1));
                zono_prs_vert = deleteAligned(project(prs_zono, 1));
                fhi = center(zono_cur_vert) +generators(zono_cur_vert);
                flo = center(zono_cur_vert) -generators(zono_cur_vert);
                phi = center(zono_prs_vert) +generators(zono_prs_vert);
                plo = center(zono_prs_vert) -generators(zono_prs_vert);
                t_now = (AH.t_move + AH.t_failsafe_move)/n*t_idx;
                sh_frs.XData = [sh_frs.XData t_now t_now];
                sh_frs.YData = [sh_frs.YData flo+agent_state(1) fhi+agent_state(1)];
                sh_prs.XData = [sh_prs.XData t_now t_now];
                sh_prs.YData = [sh_prs.YData plo+agent_state(1) phi+agent_state(1)];
            end
        end
        
        function [XData, YData] = traj_selected_zono(AH,O,K,agent_state,v_idx,ka_idx,vd_idx,theta_idx,theta_dot_idx)
            vd_idx_relative = vd_idx-v_idx+2;%convert from abs index to relative idx.
            zono_peak     = AH.zono_full.res{v_idx,ka_idx,vd_idx_relative,1};
            zono_failsafe = AH.zono_full.res{v_idx,ka_idx,vd_idx_relative,2};
            error_data_cur = AH.zono_full.error_table{v_idx,ka_idx,vd_idx_relative,theta_idx,theta_dot_idx};
            zono_one_all = [zono_peak; zono_failsafe];
            n = length(zono_one_all);
            
            %             figure(1);clf; hold on;
            pos_cart = agent_state(1);
            O_pts = O - pos_cart;
            %             yline(O_pts(1));yline(O_pts(2));yline(O_pts(4));yline(O_pts(5));
            XData = [];
            YData = [];
            %sh_prs = scatter([],[],10);
            for t_idx = 1:n
                
                zono_one = zono_one_all{t_idx}{1};
                
                if isempty(zono_one)
                    error("no such zonotope"+num2str(v_ini_idx)+num2str(y_des_idx)+num2str(h_ini_idx)+num2str(del_idx)+num2str(vd_idx));
                end
                %                 y_des = K(2);
                %                 vx_des = K(1);
                
                %                  P.vx_des,P.y_des - agent_state(2),agent_state(3),agent_state(4),agent_state(5)
                kv_soft = bound_values(agent_state(2),AH.zono_full.vi(v_idx)-AH.zono_full.kvg,AH.zono_full.vi(v_idx)+AH.zono_full.kvg);
                ka_soft = bound_values(agent_state(5),AH.zono_full.ka(ka_idx)-AH.zono_full.kag,AH.zono_full.ka(ka_idx)+AH.zono_full.kag);
                zono_one = zonotope_slice(zono_one, [3;4], [kv_soft;ka_soft]);
                %                zono_one = zonotope_slice(zono_one, [4;5], [vi_soft;agent_state(3)]);
                centererr = mean(error_data_cur(t_idx,:));
                gen_err = zeros(5,5); gen_err(1,1) = (error_data_cur(t_idx,2)- error_data_cur(t_idx,1))/2;
                err_zono  = zonotope([[centererr;0;0;0;0], gen_err]);
                
                
                full_zono = zono_one + err_zono;
                frs_zono = zonotope_slice(full_zono, 2, K);
                prs_zono = zonotope_slice(zono_one, 2, K);
                %Slice, Do it now!
                zono_cur_vert = deleteAligned(project(frs_zono, 1));
                zono_prs_vert = deleteAligned(project(prs_zono, 1));
                fhi = center(zono_cur_vert) +generators(zono_cur_vert);
                flo = center(zono_cur_vert) -generators(zono_cur_vert);
                phi = center(zono_prs_vert) +generators(zono_prs_vert);
                plo = center(zono_prs_vert) -generators(zono_prs_vert);
                t_now = (AH.t_move + AH.t_failsafe_move)/n*t_idx;
                XData = [XData t_now t_now];
                YData = [YData flo+agent_state(1) fhi+agent_state(1)];
            end
        end
        
        function [K, replaced]= find_replace_action(AH, O,k_user,agent_state,v_idx,ka_idx,vd_idx,theta_idx,theta_dot_idx)
            %select prs and ers based on initial condition
            vd_idx_relative = vd_idx-v_idx+2;%convert from abs index to relative idx.
            
            zono_peak     = AH.zono_full.res{v_idx,ka_idx,vd_idx_relative,1};
            zono_failsafe = AH.zono_full.res{v_idx,ka_idx,vd_idx_relative,2};
            error_data_cur = AH.zono_full.error_table{v_idx,ka_idx,vd_idx_relative,theta_idx,theta_dot_idx};
            zono_prs = [zono_peak; zono_failsafe];
            n = length(zono_prs);
            if n == 0
                warning("no such zonotope kv,ka,kpeak ="+num2str(v_idx)+num2str(ka_idx)+num2str(vd_idx));
                K = [];
                return
            end
            % convert obstacle to local coordinate
            pos_cart = agent_state(1);
            O_pts = O - pos_cart;
            
            % find the right place IN each prs
            % v; k_pk; kv;ka; t1
            obst_dim = [1];%integrated velocity = x
            k_dim = [2];% peak or desired velocity dimension
            A_con = {};
            b_con = {};
            s_con = {};
            buffer_dist = 0;
            
            for t_idx = 1: n
                zono_one = zono_prs{t_idx}{1};
                
                % bound the values so slicing won't have problem, max
                % values determined by simulating beforehand
                kv_soft = bound_values(agent_state(2),AH.zono_full.vi(v_idx)-AH.zono_full.kvg,AH.zono_full.vi(v_idx)+AH.zono_full.kvg);
                ka_soft = bound_values(agent_state(5),AH.zono_full.ka(ka_idx)-AH.zono_full.kag,AH.zono_full.ka(ka_idx)+AH.zono_full.kag);
                %                 kvd_soft =
                
                zono_one = zonotope_slice(zono_one, [3;4], [kv_soft;ka_soft]);
                
                % get error data from table
                %                 center_bb = error_data_cur(1:2,t_idx);
                %                 h =  error_data_cur(3,t_idx);
                %                 gen = error_data_cur(4:5,t_idx);
                %                 len = gen(1);
                %                 width = gen(2);
                %                 ego_gen = [[cos(h)*len; sin(h)*len], [sin(-h)*width; cos(-h)*width]];
                %                 gen_err = zeros(7,2);gen_err(1:2,1:2) = ego_gen;
                %                 err_zono  = zonotope([[center_bb(1);center_bb(2);0;0;0;0;0], gen_err]);
                centererr = mean(error_data_cur(t_idx,:));
                gen_err = zeros(5,5); gen_err(1,1) = (error_data_cur(t_idx,2)- error_data_cur(t_idx,1))/2;
                err_zono  = zonotope([[centererr;0;0;0;0], gen_err]);
                
                zono_one = (err_zono+zono_one);
                
                %                 z_with_err_plot = project(zono_one, [1, 2])
                %                     p4.FaceAlpha = 0.02;
                
                Z = zono_one.Z;
                A_obs_array= []; b_obs_array = [];size_array=[]; size_idx = 0;
                
                %consider each obstacle as a halfspace
                for obst_idx = 1:(size(O_pts,1)+1)
                    one_obst = O_pts((obst_idx-1)*3+1:obst_idx*3-1);
                    obst_center = mean(one_obst);
                    obst_gen = (one_obst(2)-one_obst(1))/2;
                    obst_zono  = zonotope([obst_center, obst_gen]);
                    
                    obstacle = obst_zono.Z;
                    
                    c = Z(obst_dim, 1);
                    G = Z(:, 2:end);
                    
                    for k_idx = 1:length(k_dim)
                        [~, k_col(k_idx)] = find(G(k_dim(k_idx), :) ~= 0); % find "k-sliceable" generators
                    end
                    k_slc_G = G(obst_dim, k_col);
                    k_no_slc_G = G(obst_dim, :);
                    k_no_slc_G(:, k_col) = [];
                    
                    buff_obstacle_c = [obstacle(:, 1) - c];
                    buff_obstacle_G = [obstacle(:, 2:end), k_no_slc_G, buffer_dist]; % obstacle is "buffered" by non-k-sliceable part of FRS
                    buff_obstacle_G(:, ~any(buff_obstacle_G)) = []; % delete zero columns of G
                    buff_obstacle = [buff_obstacle_c, buff_obstacle_G];
                    [A_obs, b_obs] = polytope_PH(buff_obstacle); % turn zonotope into polytope
                    A_obs_array = [A_obs_array;A_obs];
                    b_obs_array = [b_obs_array;b_obs];
                    size_idx = size_idx + length(b_obs);
                    size_array = [size_array;size_idx];
                end
                A_con{t_idx} = A_obs_array*k_slc_G; % now polytope is over coefficients of k_slc_G generators
                b_con{t_idx} = b_obs_array;
                s_con{t_idx} = size_array;
            end
            
            kpkc = AH.zono_full.kpeak(vd_idx);kpkg = AH.zono_full.kpkg;
            %             k2c = AH.y_array(y_des_idx);k2g = AH.zono_full.kyg;
            c_k = [kpkc];
            g_k = [kpkg];
            
            % grid over parameter space
            if AH.flags.discrete_flag == false
                k_sample = linspace(kpkc - kpkg, kpkc + kpkg, 9);
                k_start_idx = 5;
                %                 k2_sample = linspace(k2c - k2g, k2c + k2g, 9);
                %                 k1_usr_idx = 5; k2_usr_idx = 5; %start with the center of parameter space
            else
                k_sample = linspace(kpkc - kpkg, kpkc + kpkg, 5);
                k_start_idx = 3;
            end
            k_delta = k_sample(2)-k_sample(1);
            avaliable_action_set = zeros(0,2);
            
            for replace_idx = 0:length(k_sample)
                % was using a sprial search to make algorithm faster.
                % However, the first safe index that the sprial search
                % finds is not the closest one. So instead, just sprial
                % search the whole parameter space and record the avaliable
                % ones with distance
                if replace_idx == 0
                    K = k_user;
                else
                    kidx = k_start_idx + ceil((replace_idx - 1)/2) * (-1)^replace_idx;
                    K = k_sample(kidx);
                end
                
                %final parameter check
                if K > AH.v_max || K < -AH.v_max
                    continue;
                end
                
                
                lambdas = (K - c_k)./g_k; % given a parameter, get coefficients on k_slc_G generators
                if any(abs(lambdas) > 1)
                    continue; %When calling replace action on K_user \notin vd_idx, need this to skip the one out of range
                end
                for k = 1:length(A_con)
                    Zk_tmp = A_con{k}*lambdas - b_con{k}; % A*lambda - b <= 0 means inside unsafe set
                    safe = max(Zk_tmp(1:s_con{k}(1))) > 0;
                    for m = 1: length(s_con{k})-1
                        if ~safe
                            break %break the m loop for obstacle;
                        end
                        safe = max(Zk_tmp(s_con{k}(m)+1:s_con{k}(m+1))) > 0;
                    end
                    if ~safe
                        break % %break the k loop for time idx;
                    end
                end
                if safe
                    %                     if AH.plot_flag
                    %                         h1.XData= [h1.XData K(1)];
                    %                         h1.YData= [h1.YData K(2)];
                    %                     end
                    %save avaliable parameters in a set
                    avaliable_action_set = [avaliable_action_set; K (k_user-K)^2];
                    if replace_idx == 0
                        break;% faster if first selected action already safe
                    end
                end
                
            end
            
            if isempty(avaliable_action_set)
                K = [];
                replaced = 1;
            else
                % find the closest action based on 2nd norm
                [~,best_replace]=min(avaliable_action_set(:,2));
                K = avaliable_action_set(best_replace,1);
                if abs(k_user-K) < (k_delta - AH.eps)
                    replaced = 0;
                else
                    replaced = 1;
                end
                %                 if AH.plot_flag
                %                     figure(9);subplot(1,2,2)
                %                     %                     scatter(K(1),K(2),'g');
                %                     scatter(K(1),K(2),200,'g','filled','MarkerEdgeColor','k','LineWidth',2)
                %
                %                 end
            end
            
            
        end
        %% generate reference trajectory
        function [T,U,Z] = gen_ref(AH,action_parameterized)
            [ka,kv]=AH.get_ka_kv();
            [T,U,Z] =parameterized_cartpole_traj(action_parameterized,kv,ka);
            
            %             ylim([-8, 8])
            Z = Z+AH.A.state(1, end);
            if AH.plot_cart_flag
                figure(1);
                plot(T,U,T,Z);
            end
        end
        
        %% convert action to parameter
        function k = convert_action_to_parameter(AH,action,~)
            [~,kv] = AH.get_ka_kv();
            k = action + kv;
            if abs(k) > AH.v_max
                k = sign(k)*AH.v_max;
            end
        end
        function [ka, kv] = get_ka_kv(AH)
            agent_info = AH.get_agent_info();
            agent_state = agent_info.state(:,end);
            agent_time = agent_info.time(end);
            if size(agent_info.state,2) < 2
                previous_state = agent_state;
                previous_time = agent_time-0.001;
            else
                previous_state = agent_info.state(:,end-1);
                previous_time = agent_info.time(end-1);
            end
            ka =  (agent_state(2)- previous_state(2))/(agent_time-previous_time);
            kv     = agent_state(2);
        end
        
        %% plotting
        function plot(AH)
            %             % quadrotor_agentHelper.plot() ;
            %             %
            %             % Plot the current trajectory and corresponding FRS zonotopes.
            %
            %             % set hold on if it isn't already on
            %             hold_check = hold_switch() ;
            %
            %             % plot the current trajectory
            %             AH.plot_current_trajectory() ;
            %
            %             % plot the reachable set for the current trajectory
            %             if AH.plot_zonotope_reach_set_flag
            %                 Z = AH.current_trajectory ;
            %                 k = AH.current_parameter ;
            %
            %                 if ~isempty(Z) && ~isempty(k)
            %                     x_0 = Z(1:3,1) ;
            %                     v_0 = k(1:3) ;
            %                     a_0 = k(4:6) ;
            %                     v_peak = k(7:9) ;
            %
            %                     AH.plot_zonotope_reach_set(v_peak,x_0,v_0,a_0) ;
            %                 end
            %             end
            %
            %             % set hold off if hold wasn't on originally
            %             hold_switch(hold_check) ;
        end
        
        function plot_adjust(AH)
            %             AH.plot() ;
        end
        
        
        %% plotting helper methods
        
        %% plot trajectorygen
        function plot_current_trajectory(AH)
            %             % AH.plot_current_trajectory(t)
            %             %
            %             % Plot the current trajectory (which is in the
            %             % AH.current_trajectory property)
            %
            %             if ~isempty(AH.current_trajectory)
            %                 Z = AH.current_trajectory ;
            %
            %                 % plot current plan
            %                 if check_if_plot_is_available(AH,'trajectory')
            %                     AH.plot_data.trajectory.XData = Z(1,:) ;
            %                     AH.plot_data.trajectory.YData = Z(2,:) ;
            %                     AH.plot_data.trajectory.ZData = Z(3,:) ;
            %                 else
            %                     trajectory_data = plot3(Z(1,:),Z(2,:),Z(3,:),'b--',...
            %                         'LineWidth',AH.plot_desired_trajectory_line_width) ;
            %                     AH.plot_data.trajectory = trajectory_data ;
            %                 end
            %             end
        end
        
       function [v_peak,fval,exitflag] = trajopt_sample(AH,A_con,b_con,s_con,v_0,a_0,x_0,vd_idx,replan_start_tic)
            % create sphere at v_0
           
            % check if the quadrotor has stopped too many times and reduce
            % the number of samples if so (this can help find solutions
            % faster when stuck around too many obstacles to plan quickly)'
            % evaluate constraints
            kpkc = AH.zono_full.kpeak(vd_idx);kpkg = AH.zono_full.kpkg;
            %             k2c = AH.y_array(y_des_idx);k2g = AH.zono_full.kyg;
            lb = [AH.zono_full.kpeak(vd_idx)-AH.zono_full.kpkg];
            ub = [AH.zono_full.kpeak(vd_idx)+AH.zono_full.kpkg];
               
            k_sample = linspace(lb,ub,9);
            c_k = [kpkc];
            g_k = [kpkg];
            lambdas = (k_sample - c_k)./g_k; % given a parameter, get coefficients on k_slc_G generators
            avaliable_action_set = [];

            k_start_idx = 3;

            for replace_idx = 1:length(k_sample)

                
                kidx = k_start_idx + ceil((replace_idx - 1)/2) * (-1)^replace_idx;
                if kidx < 1
                    continue;
                end
                K = k_sample(kidx);
                
                %final parameter check
                if K > AH.v_max || K < -AH.v_max
                    continue;
                end
                
                
                lambdas = (K - c_k)./g_k; % given a parameter, get coefficients on k_slc_G generators
                if any(abs(lambdas) > 1)
                    continue; %When calling replace action on K_user \notin vd_idx, need this to skip the one out of range
                end
                for k = 1:length(A_con)
                    Zk_tmp = A_con{k}*lambdas - b_con{k}; % A*lambda - b <= 0 means inside unsafe set
                    safe = max(Zk_tmp(1:s_con{k}(1))) > 0;
                    for m = 1: length(s_con{k})-1
                        if ~safe
                            break %break the m loop for obstacle;
                        end
                        safe = max(Zk_tmp(s_con{k}(m)+1:s_con{k}(m+1))) > 0;
                    end
                    if ~safe
                        break % %break the k loop for time idx;
                    end
                end
                if safe
                    %                     if AH.plot_flag
                    %                         h1.XData= [h1.XData K(1)];
                    %                         h1.YData= [h1.YData K(2)];
                    %                     end
                    %save avaliable parameters in a set
                    avaliable_action_set = [avaliable_action_set; K];
                end
                
            end
            
            J_vals = [];
            if isempty(avaliable_action_set)
               K = [];
            else
                for i = 1:length(avaliable_action_set)
                    reward = AH.predict_cartpole_reward(avaliable_action_set(i),v_0,a_0,AH.finish_count);
                    J_vals = [J_vals; -reward];
                end
                
            end
            [~,min_idx] = min(J_vals) ;

            
            v_peak = avaliable_action_set(min_idx);
            fval = J_vals(min_idx);

            exitflag = 1 ;
       end 
        
        
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
        
        function cost = cartpole_negative_reward_cost(AH,K,kv,ka,start_tic,timeout,finishcount)
            %             dcost = zeros(2,1);
            reward= AH.predict_cartpole_reward(K,kv,ka,finishcount);
            %             agent_info = AH.get_agent_info;
            %             agent_info.state(1,end) = agent_info.state(1,end)+ x;
            %             agent_info.state(2,end) = agent_info.state(2,end)+ y;
            %             agent_info.state(4,end) = v;
            %
            %             ob = AH.S.W.get_ob(agent_info);
            %             cost = -AH.S.W.getRew(agent_info,ob);% minimizing cost, negative of rew
            cost = -reward;
            if toc(start_tic) > timeout
                AH.timeout_err_counter = AH.timeout_err_counter+1;
                error('Timed out while evaluating cost function!')
            end
        end
        
        function [c, ceq, gc, gceq] = eval_zono_cartpole_cons(AH, K, A_con, b_con, s_con,vd_idx,start_tic, timeout_t_pk)
            kpkc = AH.zono_full.kpeak(vd_idx);kpkg = AH.zono_full.kpkg;
            %             k2c = AH.y_array(y_des_idx);k2g = AH.zono_full.kyg;
            c_k = [kpkc];
            g_k = [kpkg];
            
            
            c = [];
            ceq = [];
            gc = [];
            gceq = [];
            lambdas = (K - c_k)./g_k; % given a parameter, get coefficients on k_slc_G generators
            for k = 1:length(A_con)
                Zk_tmp = A_con{k}*lambdas - b_con{k}; % A*lambda - b <= 0 means inside unsafe set
                [min_c, min_i]=min(-Zk_tmp(1:s_con{k}(1)));
                c = [c min_c]; % this and below can be combined, but lazy so leave them here
                gc= [gc -A_con{k}(min_i,:)'];
                for m = 1: length(s_con{k})-1
                    [min_c, min_i] = min(-Zk_tmp(s_con{k}(m)+1:s_con{k}(m+1)));
                    c = [c min_c];
                    gc = [gc -A_con{k}(min_i+s_con{k}(m),:)'];
                end
            end
            
%             if toc(start_tic) > timeout_t_pk
%                 AH.timeout_err_counter = AH.timeout_err_counter+1;
%                 error('Timed out while evaluating constraint function!')
%             end
        end
        function Reward = predict_cartpole_reward(AH,K,kv,ka,finishcount)
%             world_info = struct;
%             world_info.obstacles = [-5 -3 NaN 3 5];
%             stop_sim = 0;
%             start_plan_tic = tic;
            
            k = AH.convert_action_to_parameter(K,AH.flags.discrete_flag);
            %AH_copy = copy(AH);
            
            AH.stuck_count = 0;
            [T, U, Z] = AH.gen_ref(k);
            t = AH.t_move;
            z =  AH.A.state(:,end);
            
            n = 10;
            
            for i = 1:n
                state_dot = AH.A.dynamics(t/n*i,z, T, U, Z);
                z = z + t/n * state_dot;
            end
            pred_state = z;
%             Observation(1) = sin(pred_state(3));
%             Observation(2) = cos(pred_state(3));
%             Observation(3) = pred_state(4);
%             Observation(4) = pred_state(1);
%             Observation(5) = pred_state(2);
%             
%             Observation= Observation';
            
           
            theta = pred_state(3);
            
            Reward = (cos(theta)+1)/2 - 0.1*sign(pred_state(1))*sign(pred_state(2)) -0.05*abs(pred_state(4));
        end
        function [k, no_action_found] = gen_param(AH,world_info)
            agent_info = AH.get_agent_info();
            agent_state = agent_info.state(:,end); % X0;Xd0;T0;Td0
%             x_0 = agent_state(1); %X0
            O_all = world_info.obstacles; %[-5 -3 NaN 3 5]; boundary is obstacle
            no_action_found = 1;
            
            %bounds = world_info.bounds;
            O = O_all; % no heading needed to transform to local
            %% rough TRajecotry Optimization
            
            [cur_ka, kv]=AH.get_ka_kv();
            theta = agent_state(3);
            theta_dot = agent_state(4);
            
            [~,ka_idx] =min (abs(AH.zono_full.ka - cur_ka));
            %[~,vd_idx]  = min(abs(AH.zono_full.kpeak - (k_user)));%aka k peak idx
            [~,v_idx] = min(abs(AH.zono_full.vi - kv));
            theta_r = rem(theta,2*pi);%find the remainder
            if theta_r < -pi
                theta_r = theta_r + 2*pi;
            elseif theta_r > pi
                theta_r = theta_r - 2*pi;
            end
            [~,theta_idx] = min(abs([AH.zono_full.theta_a 3.1415] - theta_r));
            if theta_idx == 5
                theta_idx = 1;
            end
            %theta has circular property
            
            [~,theta_dot_idx] = min(abs(AH.zono_full.theta_dot_a - theta_dot));
            
            kv_soft = bound_values(kv,AH.zono_full.vi(v_idx)-AH.zono_full.kvg,AH.zono_full.vi(v_idx)+AH.zono_full.kvg);
            ka_soft = bound_values(cur_ka,AH.zono_full.ka(ka_idx)-AH.zono_full.kag,AH.zono_full.ka(ka_idx)+AH.zono_full.kag);
            
            agent_state = [agent_state; cur_ka];
            agent_state(2) = kv_soft;
            agent_state(5) = ka_soft;
            
            %opt_order = [1];
            %             else
            %                 opt_order = [2 1];
            %             end
            k_vec = [];
            f_valvec = [];
            timeout_t_pk = 100;
            for opt_idx = 1:3
                start_tic = tic;
                
                if opt_idx == 1
                    kpkidx = v_idx;
                elseif opt_idx == 2
                    kpkidx = v_idx - 1;
                elseif opt_idx == 3
                    kpkidx = v_idx + 1;
                end
                
                if kpkidx < 1 ||  kpkidx > length(AH.zono_full.kpeak)
                    %outside of selectable zonotopes, skip
                    continue;
                end
                
                %run optimization with cost function
                % 1. specify cost function %same as RL
                % 2. generate constraints
                [A_con,b_con,s_con] = generate_cartpole_trajopt_constraints(AH,agent_state,ka_idx,v_idx,kpkidx,theta_idx,theta_dot_idx,O);
                
                if AH.optmal_opt == "fmincon"
                    cost = @(K) AH.cartpole_negative_reward_cost(K,agent_state(2),agent_state(5),start_tic, timeout_t_pk,AH.finish_count);% one future step cost

                    cons = @(K) AH.eval_zono_cartpole_cons(K, A_con, b_con, s_con,kpkidx,start_tic, timeout_t_pk) ;% first opt_idx is for low spd, second for high spd

                    %%%%%%%
                    %                 vd_lb = AH.zono_full.v_des_range(opt_idx) -  AH.zono_full.kvdg;
                    %                 vd_ub = AH.zono_full.v_des_range(opt_idx) +  AH.zono_full.kvdg;
                    %                 vd_ub = min([vd_ub, AH.spdub]);
                    %                 y_lb  = AH.zono_full.y_range(1) - AH.zono_full.kyg;
                    %                 y_ub  = AH.zono_full.y_range(1) + AH.zono_full.kyg;

                    lb = [AH.zono_full.kpeak(kpkidx)-AH.zono_full.kpkg];
                    ub = [AH.zono_full.kpeak(kpkidx)+AH.zono_full.kpkg];
                    %                 lb =[vd_lb; y_lb];
                    %                 ub =[vd_ub; y_ub];
                    k_sample = linspace(lb,ub,10);
                    %fprintf('Current opt_idx = %d', opt_idx);
                    constraint_feasible = [];

                    for i = 1: length(k_sample)
                        constraint_feasible = [constraint_feasible any(cons(k_sample(i)) < 0)];
                    end
                    constraint_feasible;
                    initial_guess = zeros(1,1) ;



                    %                 try

                    [k,fval,exitflag,~] = fmincon(cost, initial_guess, [], [],...
                            [], [], lb, ub, cons, AH.fminconopt) ;
                elseif AH.optmal_opt == "sample"
                    [k,fval,exitflag] = trajopt_sample(AH,A_con,b_con,s_con, agent_state(2),agent_state(5),agent_state(1),kpkidx,start_tic);
                else
                    error('Wrong sampling option');
                end
                                %                 catch
                %                     exitflag = -1;%timeout
                %                 end
                if exitflag ~= 1
                    display(exitflag)
                end
                if exitflag == 1
                    k_vec = [k_vec k];
                    f_valvec = [f_valvec fval];
                end
            end
            [~,opt_idx] = min(f_valvec);
            if ~isempty(opt_idx)
                k = k_vec(:,opt_idx);
                no_action_found = 0;
            else
                display('no action found')
                k = 0;
            end
            %return by adding to class var
            %             AH.y_des = k(2)+AH.A.state(2,end);
            %             AH.vx_des = k(1);
            
        end
        
        
        
        %%
    end
end



