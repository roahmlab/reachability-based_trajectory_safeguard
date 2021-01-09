classdef zonotope_planner < planner
    properties
        %         FRS
        %         FRS_polynomial_structure
        class_num
        point_spacing
        bounds_as_obstacle
        lookahead_distance = 1.5; % default value
        current_obstacles_raw
        current_obstacles_in_FRS_coords
        y_des = 6
        y_des_old
        y_plan
        y_plan_old
        y_cur_array
        idx_iter = 1
        
        vx_des = 27
        mode = 'N'
        v_array;
        del_array;
        h_array;
        y_array;
        %zono stuff
        zono_peak;
        zono_stop;
        footprint_buffer;
        
        zono_exec = [0;0];
        zono_color = 'b';
        reset_prefer = false;
        plot_flag = 1;
        eps = 0.1;
    end
    methods
        function P = zonotope_planner(varargin)
            name = 'Zonotope planner' ;
            buffer = 0.35/2; % we'll overwrite this if we need to
            HLP = straight_line_HLP() ; % this is part of the simulator repo
            
            P = parse_args(P,'name',name,'buffer',buffer,'HLP',HLP,varargin{:}) ;
            if P.mode == 'Z'
                P.zono_peak = load('zonotpe_peak_6.12.mat');
                P.zono_stop = load('zonotpe_stop_6.12.mat');
                
                P.h_array = P.zono_peak.h_range;
                P.y_array = P.zono_peak.y_range;
                P.v_array = P.zono_peak.v_range;
                P.del_array = P.zono_peak.del_range;
                rot_max = 0.2;
                rot_len = cos(rot_max)* 4.8 + sin(rot_max)*2;
                rot_width = sin(rot_max)*4.8 + cos(rot_max)*2;
                P.footprint_buffer = zonotope([[0;0;0;0;0;0;0;0;0;0;0],diag([rot_len/2;rot_width/2;0;0;0;0;0;0;0;0;0 ])]);
            end
            
        end
        
        function setup(P,agent_info,world_info)
            
            P.y_des = agent_info.state(2)
            P. vx_des = 0;
            P.bounds = world_info.bounds;
            
            xlo = P.bounds(1) ; xhi = P.bounds(2) ;
            ylo = P.bounds(3) ; yhi = P.bounds(4) ;
            
            Blower = [xlo, xhi, xhi, xlo, xlo ; ylo, ylo, ylo-1, ylo-1, ylo] ;
            Bupper = [xlo, xhi, xhi, xlo, xlo ; yhi, yhi, yhi+1, yhi+1, yhi] ;
            B = [Blower, nan(2,1), Bupper] ;
            
            P.bounds_as_obstacle = B ;
            
            P.HLP.goal = world_info.goal ;
            P.HLP.default_lookahead_distance = P.lookahead_distance ;
            
            %             y_cur_array = ones(10, 1)*agent_info.state(2, end);
            
            P.current_plan.T = [] ;
            P.current_plan.U = [] ;
            P.current_plan.Z = [] ;
        end
        
        function [T,U,Z] = replan(P,agent_info,world_info)
            
            
            
            agent_state = agent_info.state(:,end) ; % (x,y,h,v)
            v_cur = agent_state(4) ;
            y_cur = agent_state(2) ;
            
            upper_bd = 12+6;% y
            lower_bd = 0-6;% y
            spdlb = 0;
            spdub = 10;
            
            
            if P.mode == 'R'% ||  P.mode == 'Z'
                y_des = agent_info.lane_des_action * 4 + 2;
                if abs(y_cur - y_des) > 1
                    P.y_des = y_cur + sign(y_des-y_cur);
                else
                    P.y_des = y_des;
                end
                
                P.vx_des = round((v_cur + agent_info.vx_des)/2)*2;
                if  P.vx_des > 32
                    P.vx_des= 32;
                end
                if P.vx_des < 22
                    P.vx_des = 22;
                end
            end
            
            %all modes, mostly for N     % v_cur + action [-4,2]
            P.reset_prefer = true;
            if abs(agent_info.vx_des)> P.eps
                P.vx_des = round(v_cur,1) + agent_info.vx_des;
            end
            if  P.vx_des > spdub
                P.vx_des= spdub;
            end
            if P.vx_des < spdlb
                P.vx_des = spdlb;
            end
            if abs(agent_info.y_des)> P.eps
                P.y_des = y_cur + agent_info.y_des;    % vy_cur
            end
            if  P.y_des > upper_bd
                P.y_des= upper_bd;
            end
            if P.y_des < lower_bd
                P.y_des = lower_bd;
            end
            
            if P.mode == 'Z'
                tic
                %here we should already have ready the y_des and v_des
                O_all = world_info.obstacles;
                [~,y_des_idx] =min (abs(P.y_array - (P.y_des -agent_state(2))));
                [~,v_ini_idx]  = min(abs(P.v_array - agent_state(4)));
                [~,h_ini_idx] = min(abs(P.h_array - agent_state(3)));
                [~,vd_idx] = min(abs(P.vx_des - P.v_array));
                [~,del_idx] = min(abs(P.del_array - agent_state(5)));
                %       v_array = [22 24 26 28 30 32];t
                %         h_array = deg2rad([-10: 0]);
                %         y_array = [-4:4];
                %
                %         A_MM = 0
                %         A_AM = 1
                %         A_BM = 2
                %         A_HM = 3
                %         A_ML = 4
                %         A_AL = 5% Not used at all
                %         A_BL = 6
                %         A_HL = 7 %not used often (<10)
                %         A_MR = 8
                %         A_AR = 9% not used at all
                %         A_BR = 10
                %         A_HR = 11 %not used often (<10)
                %
                %kvd = kvc + 2*(vd_idx-2)
                
                
                O = [O_all, P.bounds_as_obstacle] ;
                %v_idx y_idx h_idx vd_idx
                %check safety
                if P.plot_flag
                    rng(0)
                end
                % add world bounds as obstacle
                %                 collision_array = P.check_zono_collide_sliced(O,agent_state,v_ini_idx,y_des_idx,h_ini_idx,del_idx,vd_idx);
                collision_array = P.check_zono_collide_halfspace(O,agent_state,v_ini_idx,y_des_idx,h_ini_idx,del_idx,vd_idx);
                
                time = toc
                if false%~all( collision_array ) %unsafe
                    P.reset_prefer = true;
                    found_replacement_v = false;
                    found_replacement_y = false;
                    found_replacement_vy= false;
                    iteration_array = unique([max(P.vx_des-4,P.v_array(1)),max(P.vx_des-2,P.v_array(1)),P.vx_des,min(P.vx_des+2,P.v_array(length(P.v_array)))]);
                    for index = iteration_array
                        [~,vd_replace_idx] =  min(abs(index - P.v_array));
                        col_replace = P.check_zono_collide(O_all,agent_state,v_ini_idx,y_des_idx,h_ini_idx,vd_replace_idx);
                        if(all(col_replace))
                            found_replacement_v = true;
                            break
                        end
                    end
                    if ~found_replacement_v
                        for y_replace_idx = [5 4 6]%goes small change first
                            col_replace = P.check_zono_collide(O_all,agent_state,v_ini_idx,y_replace_idx,h_ini_idx,vd_idx);
                            if(all(col_replace))
                                found_replacement_y = true;
                                break
                            end
                        end
                    end
                    
                    if ~found_replacement_y && ~ found_replacement_v
                        col_replace_array=zeros(30,length(iteration_array),3);
                        vd_idx  = 1;
                        y_idx = 1;
                        for index = iteration_array
                            [~,vd_replace_idx] =  min(abs(index - P.v_array));
                            for y_replace_idx = [5 4 6]
                                col_replace_array(:,vd_idx,y_idx) = P.check_zono_collide(O_all,agent_state,v_ini_idx,y_replace_idx,h_ini_idx,vd_replace_idx);
                                if(all(col_replace))
                                    found_replacement_vy = true;
                                    break
                                end
                                y_idx = y_idx +1;
                            end
                            vd_idx = vd_idx +1;
                        end
                    end
                    
                    if found_replacement_v
                        delta_v = P.v_array(vd_replace_idx) - v_cur;
                        [~,idx]= min(abs(delta_v - P.delv_array));
                        P.class_num = idx-1;
                        if P.class_num < 0 || P.class_num > 4
                            error("wrong replacement action");
                        end
                        P.vx_des = P.v_array(vd_replace_idx);
                    elseif found_replacement_y
                        %                         min (abs(P.y_array - (P.y_des -agent_state(2))))
                        delta_y = P.y_array(y_replace_idx);
                        P.y_des  = delta_y+agent_state(2);
                        if delta_y == 1
                            P.class_num = 4; %ML
                        elseif delta_y == -1
                            P.class_num = 8; %MR
                        elseif delta_y == 0
                            P.class_num = 0; %MM
                        else
                            error('wrong region 2');
                        end
                    elseif found_replacement_vy
                        delta_v = P.v_array(vd_replace_idx) - v_cur;
                        [~,idx]= min(abs(delta_v - P.delv_array));
                        idx = idx-1;
                        P.vx_des = P.v_array(vd_replace_idx);
                        delta_y = P.y_array(y_replace_idx);%P.y_array(y_replace_idx)+agent_state(2) - P.y_des;
                        P.y_des  = delta_y+agent_state(2);
                        if delta_y == 1 %L
                            P.class_num = 4; %L
                        elseif delta_y == -1 %R
                            P.class_num = 8; %BR
                        elseif delta_y == 0 %M
                            P.class_num = 0; %BR
                        else
                            error('wrong region 3')
                        end
                        P.class_num = P.class_num + idx;
                        if P.class_num > 11 || P.class_num < 0
                            error('wrong region 4')
                        end
                    else
                        %                         good_idx =
                        %                         col_replace_array
                        %                         P.y_des = agent_info.lane * 4 + 2;
                        P.class_num = 0;
                        warning('no replacement act found')
                        
                    end
                    %% action 0-11
                    %do something with the replacement zono here.
                    P.zono_color = 'r';
                else
                    P.reset_prefer = true;
                    P.zono_color = 'b';
                end
                
                %                     [~,lane] = min(abs([2 6 10]-dummy_state(2)));
                %                     P.y_des = (lane -1)*4+2;
                %                     if vx_old >= 24
                %                         P.vx_des = (vx_old - 2);
                %                     else
                %                         P.vx_des = 22;
                %                     end
                
                
                %                 P.zono_exec = local_to_world(dummy_state,[[zono_cur.Vertices(:,1);zono_cur.Vertices(1,1)],[zono_cur.Vertices(:,2);zono_cur.Vertices(1,2)]]');
                
                %
                
                
                %find acc, maintain, or brake zono, check safety
                
                %                 min_score = 999999;
                %                 for con_idx=1:size(avaliable_con,1)
                %                     candidate = P.zonotope.zono_poly{avaliable_con(con_idx,1),avaliable_con(con_idx,2)};
                %                     if hit_obs(O_pts,candidate) %if safe add to avaliable
                %                         %                     candidate_con = [candidate_con; avaliable_con(con_idx,1) avaliable_con(con_idx,2)];
                %                         zono = center(P.zono_full.res{avaliable_con(con_idx,1),avaliable_con(con_idx,2)}{end}{1});
                %                         score = sum( (zono(1:2)-z_goal_local).^2);
                %                         if score < min_score
                %                             min_score = score;
                %                             best_candidate = zono;
                %                             best_can_idx = [avaliable_con(con_idx,1),avaliable_con(con_idx,2)];
                %                         end
                %                     end
                %                 end
            end
            %             if P.vx_des < 22.1
            %                 P.vx_des = 22.1;
            %             end
            
            T= [0 0.5] ;
            U =repmat( [P.vx_des;P.y_des],1, length(T));
            Z = [0 0 0 0 0;0 0 0 0 0];
            P.current_plan.T = T ;
            P.current_plan.U = U ; %use u for v_des
            P.current_plan.Z = Z; % use z for w_des
        end
        function collision_array = check_zono_collide_halfspace(P,O,agent_state,v_ini_idx,y_des_idx,h_ini_idx,del_idx,vd_idx)
            if P.plot_flag
                figure(1);
                color = rand(1,3);
            else
                color = [1 1 1];
            end
            zono_peak_one = P.zono_peak.res{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx};
            zono_stop_one = P.zono_stop.res_stop{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx};
            n = length(zono_peak_one) +length(zono_stop_one);
            collision_array = zeros(n, 1);
            
            
            P.current_obstacles = O ;
            %             [O_pts] = compute_zonotope_obs(O,...
            %                 agent_state,P.buffer,P.point_spacing) ;
            dummy_state = agent_state;
            dummy_state(3) = 0; % since the zonotope here doesn't need to be rotated.
            %             [O_pts] = compute_zonotope_obs(O,...
            %                 dummy_state,P.buffer,P.point_spacing) ;
            O_pts = world_to_local(dummy_state,O);
            % filter out points that are too far away to be reached
            %     O_FRS = crop_points_outside_region(0,0,O_FRS,1) ;
            O_pts = O_pts';
            
            if P.plot_flag
%                 for t_idx = 1: n
%                     if t_idx <= length(zono_peak_one)
%                         zono_one = zono_peak_one{t_idx}{1};
%                         color = [0 0 1];
%                     else
%                         zono_one = zono_stop_one{t_idx- length(zono_peak_one)}{1};
%                         color = [0 1 0];
%                     end
%                     zono_one = zono_one + P.footprint_buffer;
%                     if isempty(zono_one)
%                         error("no such zonotope"+num2str(v_ini_idx)+num2str(y_des_idx)+num2str(h_ini_idx)+num2str(del_idx)+num2str(vd_idx));
%                     end
%                     y_des_soft = bound_values(P.y_des - agent_state(2),-3.9, 3.9);% can't slice on edge
%                     vx_des_soft =bound_values(P.vx_des,0.1, 9.9);
%                     zono_one = zonotope_slice(zono_one, [8;9;10], [agent_state(3);agent_state(4);agent_state(5)]);
%                     
%                     figure(1)
%                     zono_cur_vert = (polygon(project(zono_one,[1 2])))';
%                     
%                     hit_obs(dummy_state, O_pts,zono_cur_vert, P.plot_flag,color);
%                     
%                 end
                for t_idx = 1: n
                    if t_idx <= length(zono_peak_one)
                        zono_one = zono_peak_one{t_idx}{1};
                        color = [0 0 1];
                    else
                        zono_one = zono_stop_one{t_idx- length(zono_peak_one)}{1};
                        color = [0 1 0];
                    end
                    zono_one = zono_one + P.footprint_buffer;
                    if isempty(zono_one)
                        error("no such zonotope"+num2str(v_ini_idx)+num2str(y_des_idx)+num2str(h_ini_idx)+num2str(del_idx)+num2str(vd_idx));
                    end
                    y_des_soft = bound_values(P.y_des - agent_state(2),-3.9, 3.9);% can't slice on edge
                    vx_des_soft =bound_values(P.vx_des,0.1, 9.9);
                    zono_one = zonotope_slice(zono_one, [8;9;10], [agent_state(3);agent_state(4);agent_state(5)]);
                    
                    figure(1)
                    zono_cur_vert_slice = (polygon(project(zonotope_slice(zono_one,[6,7],[vx_des_soft;y_des_soft]),[1 2])))';
                    hit_obs(dummy_state, O_pts,zono_cur_vert_slice, P.plot_flag,color+[0.7 0 0]);
                    
                end
                
            end
            
%             for obs_idx = 1:(size(O_pts,1)+1)/6
%                 one_obs = O_pts((obs_idx-1)*6+1:obs_idx*6-1,:);
%                 
%                 if max(one_obs(:,1)) < -50 || min(one_obs(:,1)) > 80 
%                     continue;
%                 end
%                 obs_zono = local_to_zono(one_obs);
                
%                 
                %                     for i = 1:length(Rcont)
                %                         p_FRS = plotFilled(Rcont{i}{1}, [1, 2], 'g');
                %                         p_FRS.FaceAlpha = 0.02;
                %                         p_FRS.EdgeAlpha = 0.4;
                %                     end
                
                %% create obstacle and intersect with FRS
                %                     obs_center = [0.5; 0];
                %                     obs_gen = [[0.03; 0], [0.01; 0.01]];
                %                     obs_zono = zonotope([obs_center, obs_gen]);
%                 obstacle = obs_zono.Z;
                
                % plot obstacle;
                
                % loop through zonotopes of FRS, find "k-sliceable" generators
                % and generate obstacle constraints.
                obs_dim = [1; 2]; % note that the obstacle exists in the x-y space (not theta or v)
                k_dim = [6; 7]; % note that the parameters k are in the 5th and 6th rows of the zonotopes
                buffer_dist = 0; % assume no buffer.
                
                A_con = {};
                b_con = {};
                s_con = {};
                for t_idx = 1: n
                    if t_idx <= length(zono_peak_one)
                        zono_one = zono_peak_one{t_idx}{1};
                        color = [0 0 1];
                    else
                        zono_one = zono_stop_one{t_idx- length(zono_peak_one)}{1};
                        color = [0 1 0];
                    end
                    zono_one = zono_one + P.footprint_buffer;
                    if isempty(zono_one)
                        error("no such zonotope"+num2str(v_ini_idx)+num2str(y_des_idx)+num2str(h_ini_idx)+num2str(del_idx)+num2str(vd_idx));
                    end
                    y_des_soft = bound_values(P.y_des - agent_state(2),-3.9, 3.9);% can't slice on edge
                    vx_des_soft =bound_values(P.vx_des,0.1, 9.9);
                    zono_one = zonotope_slice(zono_one, [8;9;10], [agent_state(3);agent_state(4);agent_state(5)]);
                    Z = zono_one.Z;
                    A_obs_array= []; b_obs_array = [];size_array=[]; size_idx = 0;
                    %                  P.vx_des,P.y_des - agent_state(2),agent_state(3),agent_state(4),agent_state(5)
                    for obs_idx = 1:(size(O_pts,1)+1)/6
                        one_obs = O_pts((obs_idx-1)*6+1:obs_idx*6-1,:);

                        if max(one_obs(:,1)) < -50 || min(one_obs(:,1)) > 80 
                            continue;
                        end
%                         if P.plot_flag
%                             figure(1)
%                             h=fill(O(1,(obs_idx-1)*6+1:obs_idx*6-1),O(2,(obs_idx-1)*6+1:obs_idx*6-1),'r');
%                              pause(0.2);
%                               delete(h);
%                         end
                        obs_zono = local_to_zono(one_obs);
                        obstacle = obs_zono.Z;
                    %Slice, Do it now! 6 , 7 can be sliced later, taking values vx_des_soft;y_des_soft
                    %Need to get rid of the extra zonos that have negative velocity
                       
                        c = Z(obs_dim, 1);
                        G = Z(:, 2:end);

                        [~, k_col] = find(G(k_dim, :) ~= 0); % find "k-sliceable" generators
                        k_slc_G = G(obs_dim, k_col);
                        k_no_slc_G = G(obs_dim, :);
                        k_no_slc_G(:, k_col) = [];

                        buff_obstacle_c = [obstacle(:, 1) - c];
                        buff_obstacle_G = [obstacle(:, 2:end), k_no_slc_G, buffer_dist*eye(2)]; % obstacle is "buffered" by non-k-sliceable part of FRS
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
                k1c = P.v_array(vd_idx);k1g = 0.5;
                k2c = P.y_array(y_des_idx);k2g = 2;
                c_k = [k1c; k2c];
                g_k = [k1g; k2g];
                k1_user= vx_des_soft;k2_user=y_des_soft;
                lambdas = ([k1_user; k2_user] - c_k)./g_k;
                
                
                lims = [k1c - k1g, k1c - k1g, k1c + k1g, k1c + k1g, k1c - k1g; k2c - k2g, k2c + k2g, k2c + k2g, k2c - k2g, k2c - k2g];
                    plot(lims(1, :)', lims(2, :)', 'k', 'LineWidth', 2);

                % grid over parameter space
                k1_sample = linspace(k1c - k1g, k1c + k1g, 50);
                k2_sample = linspace(k2c - k2g, k2c + k2g, 50);
                [Xk, Yk] = meshgrid(k1_sample, k2_sample);
                Zk = inf*ones(size(Xk));
                for i = 1:length(k1_sample)
                    for j = 1:length(k2_sample)
                        K = [Xk(i, j); Yk(i, j)];
                        lambdas = (K - c_k)./g_k; % given a parameter, get coefficients on k_slc_G generators
                        for k = 1:length(A_con)
                            Zk_tmp = A_con{k}*lambdas - b_con{k}; % A*lambda - b <= 0 means inside unsafe set
                            safe = true;
                            safe = max(Zk_tmp(1:s_con{k}(1))) > 0;
                            for m = 1: length(s_con{k})-1
                                if ~safe
                                    break
                                end
                              
                                safe = max(Zk_tmp(s_con{k}(m)+1:s_con{k}(m+1))) > 0;
                            end
                            if safe
                                Zk(i, j) = min(Zk(i, j), 1); % take smallest max. if it's <=0, then unsafe
                            else
                                Zk(i, j) =-1;
                            end
                        end
                    end
                end
                
                if P.plot_flag
                    % plot parameters and unsafe set
                    figure(9); clf; hold on;
                    title('Set of parameters (red is unsafe)', 'FontSize', 24);
                    xlabel('$k_1$', 'Interpreter', 'latex', 'FontSize', 24);
                    ylabel('$k_2$', 'Interpreter', 'latex', 'FontSize', 24);

                    p_unsafe = contourf(Xk, Yk, -Zk, [0, 0], 'FaceColor', 'r'); % show zero level set contours
                    
                    %% have user select a point, see the corresponding slice of FRS:
                    
                    k1_user= vx_des_soft;k2_user=y_des_soft;
                    
                    plot(k1_user, k2_user, 'b.', 'MarkerSize', 30, 'LineWidth', 6);
                    figure(1)
                    
                    
                    pause(0.1);
%                     delete(h);
                end
                
                
%             end
            
            
        end
        function collision_array = check_zono_collide_sliced(P,O,agent_state,v_ini_idx,y_des_idx,h_ini_idx,del_idx,vd_idx)
            
            % following feature for zonotope
            if P.plot_flag
                figure(1);
                color = rand(1,3);
            else
                color = [1 1 1];
            end
            zono_peak_one = P.zono_peak.res{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx};
            zono_stop_one = P.zono_stop.res_stop{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx};
            n = length(zono_peak_one) +length(zono_stop_one);
            collision_array = zeros(n, 1);
            for t_idx = 1: n
                if t_idx <= length(zono_peak_one)
                    zono_one = zono_peak_one{t_idx}{1};
                    color = [0 0 1];
                else
                    zono_one = zono_stop_one{t_idx- length(zono_peak_one)}{1};
                    color = [0 1 0];
                end
                if isempty(zono_one)
                    error("no such zonotope"+num2str(v_ini_idx)+num2str(y_des_idx)+num2str(h_ini_idx)+num2str(del_idx)+num2str(vd_idx));
                end
                %                  P.vx_des,P.y_des - agent_state(2),agent_state(3),agent_state(4),agent_state(5)
                y_des_soft = bound_values(P.y_des - agent_state(2),-3.9, 3.9);% can't slice on edge
                vx_des_soft =bound_values(P.vx_des,0.1, 9.9);
                zono_one = zonotope_slice(zono_one, [6;7;8;9;10], [vx_des_soft;y_des_soft;agent_state(3);agent_state(4);agent_state(5)]);
                %Slice, Do it now!
                zono_cur_vert = (polygon(project(zono_one+P.footprint_buffer,[1 2])))';
                %Need to get rid of the extra zonos that have negative velocity
                % zono_cur_vert = [[zono_cur.Vertices(:,1);zono_cur.Vertices(1,1)],[zono_cur.Vertices(:,2);zono_cur.Vertices(1,2)]];
                
                % buffer and discretize obstacles
                %             O_pts = buffer_polygon_obstacles(O,P.buffer,2) ;
                %                     O_pts = buffer_polygon_obstacles(O,0.1,2) ;% don't buffer since zono already buffered
                P.current_obstacles = O ;
                %             [O_pts] = compute_zonotope_obs(O,...
                %                 agent_state,P.buffer,P.point_spacing) ;
                dummy_state = agent_state;
                dummy_state(3) = 0; % since the zonotope here doesn't need to be rotated.
                [O_pts] = compute_zonotope_obs(O,...
                    dummy_state,P.buffer,P.point_spacing) ;
                
                
                % save obstacles
                %             P.current_obstacles_raw = O ; % input from the world
                % buffered and discretized
                
                %DO NOTHING
                
                % choose the current state corresponding zonotope and the
                % current desired v and desired y zonotope, see if it is
                % safe, if not, check maintain or break and do one of
                % those.
                
                %find the zono corresponding to current state/ desired
                %location.
                %                     if t_idx <= 20
                %                         local_plot_flag = true;
                %                     else
                %                         local_plot_flag = false;
                %                     end
                collision_array(t_idx) = hit_obs(dummy_state, O_pts,zono_cur_vert, P.plot_flag,color);
                %                 pause(0.05)
            end
        end
        function collision_array = check_zono_collide(P,O_all,agent_state,v_ini_idx,y_des_idx,h_ini_idx,vd_idx)
            collision_array = zeros(30, 1);
            % following feature for zonotope
            if P.plot_flag
                figure(1);
                color = rand(1,3);
            else
                color = [1 1 1];
            end
            for t_idx = 1: 30
                
                zono_cur = P.zonotope.zono_poly{v_ini_idx ,y_des_idx,h_ini_idx,vd_idx,t_idx};
                zono_cur_vert = [[zono_cur.Vertices(:,1);zono_cur.Vertices(1,1)],[zono_cur.Vertices(:,2);zono_cur.Vertices(1,2)]];
                O = [O_all{t_idx}, nan(2,1), P.bounds_as_obstacle] ;
                
                % buffer and discretize obstacles
                %             O_pts = buffer_polygon_obstacles(O,P.buffer,2) ;
                %                     O_pts = buffer_polygon_obstacles(O,0.1,2) ;% don't buffer since zono already buffered
                P.current_obstacles = O ;
                %             [O_pts] = compute_zonotope_obs(O,...
                %                 agent_state,P.buffer,P.point_spacing) ;
                dummy_state = agent_state;
                dummy_state(3) = 0; % since the zonotope here doesn't need to be rotated.
                [O_pts] = compute_zonotope_obs(O,...
                    dummy_state,P.buffer,P.point_spacing) ;
                
                
                % save obstacles
                %             P.current_obstacles_raw = O ; % input from the world
                % buffered and discretized
                
                %DO NOTHING
                
                % choose the current state corresponding zonotope and the
                % current desired v and desired y zonotope, see if it is
                % safe, if not, check maintain or break and do one of
                % those.
                
                %find the zono corresponding to current state/ desired
                %location.
                %                     if t_idx <= 20
                %                         local_plot_flag = true;
                %                     else
                %                         local_plot_flag = false;
                %                     end
                local_plot_flag = false;
                collision_array(t_idx) = hit_obs(dummy_state, O_pts,zono_cur_vert, local_plot_flag,color);
                %                 pause(0.05)
            end
            
        end
    end
end
