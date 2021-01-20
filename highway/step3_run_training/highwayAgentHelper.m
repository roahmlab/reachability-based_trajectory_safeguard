classdef highwayAgentHelper < agentHelper
    %% This class inherit the agentHelper class to help adjusting parameter choices RL makes, contains agent as a class variable
    properties
        % hard reference bounds for parameters
        upper_bd = 10;% y
        lower_bd = 2;% y
        spdlb = 1;% vx
        spdub = 4;%vx
        
        % epsilon for boundary check
        eps = 0.001;
        
        
        % reference for a single sampling time, needs this so there is a
        % old reference when a new reference comes in
        vx_des = 1; % default to 1m/s
        y_des;
        
        % from FRS, to determine initial condition and desired condition
        h_array
        y_array
        v_array
        del_array
        
        plot_flag
        %get rid of this eventually
        
        HLP = [];% high level planner, is empty for the Z and N mode, R mode gets constructed in constructor
        fminconopt = [];
        
        M; %temp dictionary storing FRS with better format
        
        % reference data for plot
        ref_Z = [];
        FRS_plotting_param = [];
        proposed_ref_Z = [];
        t_real_start = [];
        t_proposed_start = [];
        
%         timeout_err_counter = 0;
        S = [];%simulator pointer so eaiser to access world data
    end
    %% methods
    methods
        function AH = highwayAgentHelper(A,FRS_path,HLP,varargin)
            AH@agentHelper(A,FRS_path,varargin{:});
            AH.y_des = AH.A.state(2,end);
            AH.h_array = AH.zono_full.h_range; % array of range of acceptable parameter values
            AH.y_array = AH.zono_full.y_range;
            AH.v_array = AH.zono_full.v_range;
            AH.zono_full.v_range = AH.zono_full.v_range;
            AH.del_array = AH.zono_full.del_range;
            AH.HLP = HLP;
%             if ~isempty(HLP)
                o = optimoptions('fmincon') ;
                o.OptimalityTolerance = 1e-3 ;
                o.MaxIterations = 100000 ;
                o.MaxFunctionEvaluations = 10000 ;
                o.SpecifyConstraintGradient = true ;
                o.SpecifyObjectiveGradient = isa(AH.HLP,'highway_HLP');
                o.Display = 'off';
                o.CheckGradients = false;
%             end
            AH.fminconopt = o;
            % zono_full and M are essentially the same data, M is in a cleaner format, outputed in step 1
            AH.M = load("zono_full_7.13_1spd_cleaned.mat");
            AH.M = AH.M.M;
        end
        
        
        
        %% methods called in parent class
        function [K, replace_distance, replaced_flag] = adjust(AH,K_user,world_info)
            % takes a user input parameter and obstacles in world info,
            % find a nearest action that is safe.(may be itself), also
            % return the distance in parameter space and return a replace
            % flag based on if a proper replace action is found
            AH.gen_ref(K_user,0);
            replace_distance = 0;
            agent_info = AH.get_agent_info();
            agent_state = agent_info.state(:,end);
            O_all = world_info.obstacles;
            [~,y_des_idx] =min (abs(AH.y_array - (K_user(2))));
            [~,v_ini_idx]  = min(abs(AH.v_array - agent_state(4)));
            [~,h_ini_idx] = min(abs(AH.h_array - agent_state(3)));
            [~,vd_idx] = min(abs(K_user(1) - AH.v_array));
            [~,del_idx] = min(abs(AH.del_array - agent_state(5)));
            % find the index correesponding to the initial condition (v_ini, h_ini, delta_ini)
            % paramter range and desired parameter range(y_desired, and vd).
            
            % inefficient to generate bounds obstacle everytime
            bounds = world_info.bounds;
            
            xlo = bounds(1) ; xhi = bounds(2) ;
            ylo = bounds(3) ; yhi = bounds(4) ;
            
            Blower = [xlo, xhi, xhi, xlo, xlo ; ylo, ylo, ylo-1, ylo-1, ylo] ;
            Bupper = [xlo, xhi, xhi, xlo, xlo ; yhi, yhi, yhi+1, yhi+1, yhi] ;
            B = [Blower, nan(2,1), Bupper] ;% make top and bottom bounds into obstacles
            
            O = [O_all, B] ;
            
            %find the action to replace starting from current zonotope
            %index
            [K,replaced_flag] = AH.find_replace_action_global(O,agent_state,v_ini_idx,y_des_idx,h_ini_idx,del_idx,vd_idx);
            if AH.plot_flag
                if ~isempty(K)
                    AH.plot_zono_collide_sliced(O,agent_state,v_ini_idx,y_des_idx,h_ini_idx,del_idx,vd_idx,K);
                    plotting_param = struct;
                    plotting_param.O = O;
                    plotting_param.agent_state = agent_state;
                    plotting_param.v_ini_idx = v_ini_idx;
                    plotting_param.y_des_idx = y_des_idx;
                    plotting_param.h_ini_idx = h_ini_idx;
                    plotting_param.del_idx = del_idx;
                    plotting_param.vx_idx = vd_idx;
                    AH.FRS_plotting_param= [AH.FRS_plotting_param; plotting_param];
                end
            end
            if isempty(K)
                replaced_flag = 2;
            end
            if replaced_flag == 1
                replace_distance = sum((K_user - K).^2);
            end
            
        end
        
        %%
        function [k, no_action_found] = gen_param(AH,world_info)
            %return parameter based on optimization result
            % since optimization based options can sometimes take very long, error out and catch that error and return no action found
            agent_info = AH.get_agent_info();
            agent_state = agent_info.state(:,end);
            x_0 = agent_state(1:2);
            O_all = world_info.obstacles;
            
            no_action_found = 1;
            % find the index correesponding to the initial condition (v_ini, h_ini, delta_ini)
            % paramter range and desired parameter range(y_desired, and vd).
            
            % inefficient to generate bounds obstacle everytime
            bounds = world_info.bounds;
            
            xlo = bounds(1) ; xhi = bounds(2) ;
            ylo = bounds(3) ; yhi = bounds(4) ;
            
            Blower = [xlo, xhi, xhi, xlo, xlo ; ylo, ylo, ylo-1, ylo-1, ylo] ;
            Bupper = [xlo, xhi, xhi, xlo, xlo ; yhi, yhi, yhi+1, yhi+1, yhi] ;
            B = [Blower, nan(2,1), Bupper] ;% make top and bottom bounds into obstacles
            
            O = [O_all, B] ;
            dummy_state = agent_state;
            dummy_state(3) = 0;
            O = world_to_local(dummy_state,O)';
            %find the action to replace starting from current zonotope
            %index
            %% rough TRajecotry Optimization
            lookahead_distance = AH.t_move*AH.spdub;
            if isa(AH.HLP,'highway_HLP')
                d = zeros(3,1);
                [d(3),~,~,~] = AH.S.W.getCar_abs(2,agent_info);%top lane
                [d(2),~,~,~] = AH.S.W.getCar_abs(4,agent_info);% middle lane
                [d(1),~,~,~] = AH.S.W.getCar_abs(6,agent_info);%bottom lane
                
                x_des = AH.HLP.get_waypoint(agent_info,world_info,d) ;
                % local coordinate desired
                x_des = x_des - x_0 ;
            else
                %direct -rew optimization, no HLP and no x_des
            end
            %find all options
            [~,v_ini_idx]  = min(abs(AH.v_array - agent_state(4)));
            [~,h_ini_idx] = min(abs(AH.h_array - agent_state(3)));
            [~,del_idx] = min(abs(AH.del_array - agent_state(5)));
            %y_idx = 1;
            % desired vx idx is embedded in the table we can choose using
            % opt_idx
            u0  = AH.v_array(v_ini_idx);
            hi  = AH.h_array(h_ini_idx);
            del = AH.del_array(del_idx);
            %             tb  = AH.M(char("u0="+num2str(u0)+"hi="+num2str(hi)+"delta="+num2str(del)+"_tb"));
            zono_ini = AH.M(char("u0="+num2str(u0)+"hi="+num2str(hi)+"delta="+num2str(del)+"_zono"));
            %             d = vecnorm(repmat(x_des,1,2) - tb(2:3,:));
            % Make sure states are in bounds
            vi_soft =bound_values(agent_state(4),AH.zono_full.v_range(v_ini_idx)-AH.zono_full.kvig+AH.eps, AH.zono_full.v_range(v_ini_idx)+AH.zono_full.kvig-AH.eps);
            h_soft = bound_values(agent_state(3),AH.zono_full.h_range(h_ini_idx)-AH.zono_full.h_ini_range+AH.eps,AH.zono_full.h_range(h_ini_idx)+AH.zono_full.h_ini_range-AH.eps);
            agent_state(4) = vi_soft;
            agent_state(3) = h_soft;
            
            %             [~,spd_idx]=min(d);
            %             if spd_idx == 1
            opt_order = [1 2];
            %             else
            %                 opt_order = [2 1];
            %             end
            k_vec = [];
            f_valvec = [];
            timeout_t_pk = 2;
            for opt_idx = opt_order %time them seperately since it can be run in parallel
                start_tic = tic;
                %run optimization with cost function
                % 1. specify cost function %same as RL
                % 2. generate constraints
                [A_con,b_con,s_con] = generate_highway_trajopt_constraints(agent_state,zono_ini{opt_idx}, O);
                % 3. Fmincon with constraints
                % 3.1 define cost fun (the 2 second future distance to the waypoint)
                % cost_type 1: predicted position to waypoint
                if isa(AH.HLP,'highway_HLP')
                    cost = @(K) highway_cost(K, agent_state(4),agent_state(3),x_des,start_tic, timeout_t_pk) ;% one future step cost
                else %directly optimize for reward 
                    cost = @(K) AH.highway_negative_reward_cost(K,agent_state(4),agent_state(3),start_tic, timeout_t_pk);% one future step cost
                    % cost_type 2: -RL reward shortsighted
                end
                
                
                cons = @(K) AH.eval_zono_highway_cons(K, A_con, b_con, s_con,opt_idx,start_tic, timeout_t_pk) ;% first opt_idx is for low spd, second for high spd
                %                 cons = @(v_peak)  ;
                
                vd_lb = AH.zono_full.v_des_range(opt_idx) -  AH.zono_full.kvdg;
                vd_ub = AH.zono_full.v_des_range(opt_idx) +  AH.zono_full.kvdg;
                vd_ub = min([vd_ub, AH.spdub]);
                y_lb  = AH.zono_full.y_range(1) - AH.zono_full.kyg;
                y_ub  = AH.zono_full.y_range(1) + AH.zono_full.kyg;
                
                lb =[vd_lb; y_lb];
                ub =[vd_ub; y_ub];
                
                initial_guess = zeros(2,1) ;
                
                try
                    [k,fval,exitflag,~] = fmincon(cost, initial_guess, [], [],...
                        [], [], lb, ub, cons,AH.fminconopt) ;
                catch
                    exitflag = -1;%timeout
                end
                
                if exitflag == 1
                    k_vec = [k_vec k];
                    f_valvec = [f_valvec fval];
                end
            end
            % find the min out of the two avaliable desired velocity range.
            [~,opt_idx] = min(f_valvec);
            if ~isempty(opt_idx)
                k = k_vec(:,opt_idx);
                no_action_found = 0;
            else
                k = [0;0];
            end
            %return by adding to class var
            AH.y_des = k(2)+AH.A.state(2,end); % since y is in local coordinate, change to global coordinate.
            AH.vx_des = k(1);
            
        end
        
        %%
        function [cost ] = highway_negative_reward_cost(AH,K,kv,kh,start_tic,timeout)
%             dcost = zeros(2,1);
            [x,y,v] = predict_highway_position( K,kv,kh );
            agent_info = AH.get_agent_info;
            agent_info.state(1,end) = agent_info.state(1,end)+ x;
            agent_info.state(2,end) = agent_info.state(2,end)+ y;
            agent_info.state(4,end) = v;
            
            ob = AH.S.W.get_ob(agent_info);
            cost = -AH.S.W.getRew(agent_info,ob);% minimizing cost, negative of rew
            if toc(start_tic) > timeout
                AH.timeout_err_counter = AH.timeout_err_counter+1;
                error('Timed out while evaluating cost function!')
            end
        end
        %%
        function K = convert_action_to_parameter(AH,action,discrete_flag)
            % called in parent class. differnet for each helper
            delta_struct = struct;
                delta_struct.vx_des = (action(1)+1)/2*6-4;   %%%%%%  vx
                delta_struct.y_des = action(2)*1;%+1)/2*8+2;    %%%%%%  y
%             end
            
            if delta_struct.vx_des > 2
                delta_struct.vx_des = 2;
            elseif delta_struct.vx_des < -4
                delta_struct.vx_des = -4;
            end
            
            if delta_struct.y_des > 1
                delta_struct.y_des = 1;
            elseif delta_struct.y_des < -1
                delta_struct.y_des = -1;
            end
            
            %above is for delta vx and delta y, the following is for actual
            %parameters,      : vx and delta y
            K = AH.update_desired_parameters(delta_struct);
            
        end
        function [T, U, Z]=gen_ref(AH, K, real_reference_flag)
            % generate reference based on parameter and states
            agent_info = AH.get_agent_info();
            agent_state = agent_info.state(:,end);
            v_cur = agent_state(4) ;
            h_cur = agent_state(3) ;
            y_cur = agent_state(2) ;
            x_cur = agent_state(1) ;
            [T, U,Z] = parameterized_traj_1d_with_brake(K(1), v_cur,h_cur,K(2));
            U(2,:) = U(2,:)+y_cur;% local coordinate to global
            if ~exist('real_reference_flag','var')
                real_reference_flag = 1;
            end
            if real_reference_flag
                AH.ref_Z=[AH.ref_Z;x_cur+Z(1,:);y_cur+Z(2,:)];% for plotting
                AH.t_real_start = [AH.t_real_start;AH.A.time(end)];
            else
                AH.proposed_ref_Z=[AH.proposed_ref_Z;x_cur+Z(1,:);y_cur+Z(2,:)];% for plotting
                AH.t_proposed_start = [AH.t_proposed_start;AH.A.time(end)];
            end
            
            
        end
        
        %% helper functions
        function parameters=update_desired_parameters(AH,delta_struct)
            parameters = zeros(2,1);
            agent_info = AH.get_agent_info();
            agent_state = agent_info.state(:,end);
            v_cur = agent_state(4) ;
            %             h_cur = agent_state(3) ;
            y_cur = agent_state(2) ;
            %             x_cur = agent_state(1) ;
            %             replaced_flag = 0;
            
            if abs(delta_struct.vx_des)> AH.eps
                AH.vx_des = round(v_cur,1) + delta_struct.vx_des;
            end
            if  AH.vx_des > AH.spdub
                AH.vx_des= AH.spdub;
            end
            if AH.vx_des < AH.spdlb
                AH.vx_des = AH.spdlb;
            end
            if abs(delta_struct.y_des)> AH.eps
                AH.y_des = y_cur + delta_struct.y_des;    % vy_cur
            end
            if  AH.y_des > AH.upper_bd
                AH.y_des= AH.upper_bd;
            end
            if AH.y_des < AH.lower_bd
                AH.y_des = AH.lower_bd;
            end
            parameters(1) = AH.vx_des;
            parameters(2) = AH.y_des-y_cur;
        end
        function plot(AH)
            hold_check = false ;
            if ~ishold
                hold_check = true ;
                hold on
            end
            if ~isempty(AH.HLP)
                AH.HLP.plot;
            end
            if ~isempty(AH.proposed_ref_Z)
                plot(AH.proposed_ref_Z(end-1,:),AH.proposed_ref_Z(end,:),'k-','LineWidth',3);
                plot(AH.proposed_ref_Z(end-1,:),AH.proposed_ref_Z(end,:),'Color','y','LineWidth',3,'LineStyle','--');
                
                %                 xlim([AH.ref_Z(1,1)-20,AH.ref_Z(1,1)+30]);
            end
            if ~isempty(AH.ref_Z)
                plot(AH.ref_Z(end-1,:),AH.ref_Z(end,:),'Color',[0 0 0],'LineStyle','-','LineWidth',3);
                %                 plot(AH.ref_Z(end-1,:),AH.ref_Z(end,:),'g--','LineWidth',2);
                
                %                 xlim([AH.ref_Z(1,1)-20,AH.ref_Z(1,1)+30]);
            end
            
            if hold_check
                hold off ;
            end
            
        end
        
        function plot_adjust(AH)
            error('adjust plot not implemented')
        end
        
        function reset(AH,flags,eps_seed)
            if ~exist('eps_seed','var')
                AH.A.reset();
            else
                rng(eps_seed)
                AH.A.reset();
            end
            if isa(AH.HLP,'RRT_star_HLP')
                AH.HLP.reset_tree;
            end
            AH.y_des = AH.A.state(2,end);
            AH.vx_des = AH.A.state(4,end);
            AH.flags = flags;
            AH.ref_Z = [];
            AH.proposed_ref_Z = [];
            AH.t_real_start = [];
            AH.t_proposed_start = [];
        end
        function plot_zono_collide_sliced(AH,O,agent_state,v_ini_idx,y_des_idx,h_ini_idx,del_idx,vd_idx,K)
            % This was the old method to check collision, now used for plotting and debugging
            if AH.plot_flag
                color = [0 1 0];
            else
                color = [1 1 1];
            end
            zono_peak_peak = AH.zono_full.res{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx,1};
            zono_peak_mid = AH.zono_full.res{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx,2};
            zono_peak_stop = AH.zono_full.res{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx,3};
            error_data_cur = AH.zono_full.error_table{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx};
            zono_one_all = [zono_peak_peak; zono_peak_mid;zono_peak_stop];
            n = length(zono_one_all);
            collision_array = zeros(n, 1);
            for t_idx = 1:12:n
                
                zono_one = zono_one_all{t_idx}{1};
                
                if isempty(zono_one)
                    error("no such zonotope"+num2str(v_ini_idx)+num2str(y_des_idx)+num2str(h_ini_idx)+num2str(del_idx)+num2str(vd_idx));
                end
                y_des = K(2);
                vx_des = K(1);
                
                %                  P.vx_des,P.y_des - agent_state(2),agent_state(3),agent_state(4),agent_state(5)
                y_des_soft = bound_values(y_des,AH.zono_full.y_range(y_des_idx)-AH.zono_full.kyg+AH.eps, AH.zono_full.y_range(y_des_idx)+AH.zono_full.kyg-AH.eps);% can't slice on edge
                vx_des_soft =bound_values(vx_des,AH.zono_full.v_range(vd_idx)-AH.zono_full.kvdg+AH.eps, AH.zono_full.v_range(vd_idx)+AH.zono_full.kvdg-AH.eps);
                vi_soft =bound_values(agent_state(4),AH.zono_full.v_range(v_ini_idx)-AH.zono_full.kvig+AH.eps, AH.zono_full.v_range(v_ini_idx)+AH.zono_full.kvig-AH.eps);
                %                 sliced_zono = zonotope_slice(zono_one, [3;4;5;6], [vx_des_soft;vi_soft;agent_state(3);y_des_soft]);
                heading_soft = bound_values(agent_state(3), AH.zono_full.h_range(h_ini_idx)-AH.zono_full.h_ini_range+AH.eps, AH.zono_full.h_range(h_ini_idx)+AH.zono_full.h_ini_range-AH.eps);
                zono_one = zonotope_slice(zono_one, [4;5], [vi_soft;heading_soft]);
                %                zono_one = zonotope_slice(zono_one, [4;5], [vi_soft;agent_state(3)]);
                
                center_bb = error_data_cur(1:2,t_idx);
                h =  error_data_cur(3,t_idx);
                gen = error_data_cur(4:5,t_idx);
                len = gen(1);
                width = gen(2);
                ego_gen = [[cos(h)*len; sin(h)*len], [sin(-h)*width; cos(-h)*width]];
                gen = zeros(7,2);gen(1:2,1:2) = ego_gen;
                err_zono  = zonotope([[center_bb(1);center_bb(2);0;0;0;0;0], gen]);
                
                full_zono = zono_one + err_zono;
                sliced_zono = zonotope_slice(full_zono, [3;6], [vx_des_soft;y_des_soft]);
                %Slice, Do it now!
                zono_cur_vert = polygon(project(sliced_zono, [1, 2]))';

                dummy_state = agent_state;
                dummy_state(3) = 0; % since the zonotope here doesn't need to be rotated.
                O_pts = world_to_local(dummy_state,O);
                % filter out points that are too far away to be reached
                %     O_FRS = crop_points_outside_region(0,0,O_FRS,1) ;
                O_pts = O_pts';
                
                plot_flag = AH.plot_flag;
                %                 plot_flag = 0;
                collision_array(t_idx) = hit_obs(dummy_state, O_pts,zono_cur_vert, plot_flag,color);
            end
        end
        
        function [K, replaced] = find_replace_action_global(AH,O,agent_state,v_ini_idx,y_des_idx,h_ini_idx,del_idx,vd_idx)
            %for each combination of v_des_index and y_des_idx, try finding
            %the replacment action. start with the current index, if cannot
            %find coninue onto near ones.
            % Interestingly, there is only one combination now. v = [1 3]
            %y = [-1 1], v =[3 5] was disabled due to a bug, but that makes
            %everything simpler so the method is more clear.but maybe not
            %to its full potential.
            for replace_idx = 1:2*2*length(AH.zono_full.v_range)*length(AH.zono_full.y_range)% 2 times making sure all zonotopes are looked at
                [d1,d2] = sprial_index(replace_idx); % use a spiral indexing way to check nearby zonotopes
                k1idx = vd_idx + d1;
                k2idx = y_des_idx + d2;
                if k1idx > length(AH.zono_full.v_range)  || k1idx < 1 || k2idx < 1 ||  k2idx > length(AH.zono_full.y_range)
                    %outside of selectable zonotopes, skip
                    continue;
                end
                if AH.plot_flag
                    AH.check_zono_collide_halfspace(O,agent_state,v_ini_idx,k2idx,h_ini_idx,del_idx,k1idx);% check for full space safety, for good visualization
                end
                [K, inner_replaced_flag] = AH.find_replace_action(O,agent_state,v_ini_idx,k2idx,h_ini_idx,del_idx,k1idx);% find a action that works
                if ~isempty(K)
                    %P.check_zono_collide_sliced(O,agent_state,v_ini_idx,k2idx,h_ini_idx,del_idx,k1idx,K);
                    break
                end
            end
            if replace_idx == 1 && inner_replaced_flag == 0
                replaced = 0;
            else
                replaced = 1;
            end
        end
        function [c, ceq, gc, gceq] = eval_zono_highway_cons(AH,K,A_con,b_con,s_con,vd_idx,start_tic,timeout)
            k1c = AH.v_array(vd_idx);k1g = AH.zono_full.kvdg;
            k2c = AH.y_array(1);k2g = AH.zono_full.kyg;% y_des only 1 option
            c_k = [k1c; k2c];
            g_k = [k1g; k2g];
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
            if toc(start_tic) > timeout
                AH.timeout_err_counter = AH.timeout_err_counter+1;
                error('Timed out while evaluating constraint function!')
            end
            %we need max(Ax -b)> 0, since fmincon requires nonlin<=0
            %we specify constraint as min(b-Ax)<=0
        end
        function [K, replaced] = find_replace_action(AH,O,agent_state,v_ini_idx,y_des_idx,h_ini_idx,del_idx,vd_idx)
            % replace action based on obstacle check. takes intial
            % condition index and desired condition index, output the
            % replaced parameter and if it is replaced. return [[],true]
            % if no replacement is found.
            zono_peak_peak = AH.zono_full.res{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx,1};
            zono_peak_mid = AH.zono_full.res{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx,2};
            zono_peak_stop = AH.zono_full.res{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx,3};
            error_data_cur = AH.zono_full.error_table{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx};
            
            zono_peak_mid_stop = [zono_peak_peak; zono_peak_mid;zono_peak_stop];
            n = length(zono_peak_mid_stop);
            if n == 0
                warning("no such zonotope"+num2str(v_ini_idx)+num2str(y_des_idx)+num2str(h_ini_idx)+num2str(del_idx)+num2str(vd_idx));
                K = [];
                return
            end
            
            dummy_state = agent_state;
            dummy_state(3) = 0; % since the zonotope here doesn't need to be rotated.
            
            O_pts = world_to_local(dummy_state,O);
            
            O_pts = O_pts';
            
            
            obs_dim = [1; 2]; % note that the obstacle exists in the x-y space (not theta or v)
            k_dim = [3; 6]; % note that the parameters k are in the 3th and 6th rows of the zonotopes
            buffer_dist = 0; % assume no buffer.
            %
            A_con = {};
            b_con = {};
            s_con = {};
            for t_idx = 1: n
                
                zono_one = zono_peak_mid_stop{t_idx}{1};
                
                % bound the values so slicing won't have problem, max
                % values determined by simulating beforehand
                y_des_soft = bound_values(AH.y_des - agent_state(2),AH.zono_full.y_range(y_des_idx)-AH.zono_full.kyg+AH.eps, AH.zono_full.y_range(y_des_idx)+AH.zono_full.kyg-AH.eps);% can't slice on edge
                vx_des_soft =bound_values(AH.vx_des,AH.zono_full.v_range(vd_idx)-AH.zono_full.kvdg+AH.eps, AH.zono_full.v_range(vd_idx)+AH.zono_full.kvdg-AH.eps);
                vi_soft =bound_values(agent_state(4),AH.zono_full.v_range(v_ini_idx)-AH.zono_full.kvig+AH.eps, AH.zono_full.v_range(v_ini_idx)+AH.zono_full.kvig-AH.eps);
                h_soft = bound_values(agent_state(3),AH.zono_full.h_range(h_ini_idx)-AH.zono_full.h_ini_range+0.001,AH.zono_full.h_range(h_ini_idx)+AH.zono_full.h_ini_range-0.001);
                zono_one = zonotope_slice(zono_one, [4;5], [vi_soft;h_soft]);
                
                % get error data from table
                center_bb = error_data_cur(1:2,t_idx);
                h =  error_data_cur(3,t_idx);
                gen = error_data_cur(4:5,t_idx);
                len = gen(1);
                width = gen(2);
                ego_gen = [[cos(h)*len; sin(h)*len], [sin(-h)*width; cos(-h)*width]];
                gen_err = zeros(7,2);gen_err(1:2,1:2) = ego_gen;
                err_zono  = zonotope([[center_bb(1);center_bb(2);0;0;0;0;0], gen_err]);
                
                
                zono_one = (err_zono+zono_one);
                
                
                Z = zono_one.Z;
                A_obs_array= []; b_obs_array = [];size_array=[]; size_idx = 0;
                
                %consider each obstacle as a halfspace
                for obs_idx = 1:(size(O_pts,1)+1)/6
                    one_obs = O_pts((obs_idx-1)*6+1:obs_idx*6-1,:);
                    
                    
                    if max(one_obs(:,1)) < -20 || min(one_obs(:,1)) > 50
                        %                     if(max(one_obs(:,1)) - min(one_obs(:,1)) ) < 100 || max(one_obs(:,2))>0
                        
                        continue;
                    end
                    %                         if P.plot_flag
                    %                             figure(1)
                    %                             h=(O(1,(obs_idx-1)*6+1:obs_idx*6-1),O(2,(obs_idx-1)*6+1:obs_idx*6-1),'r');
                    %                              pause(0.2);
                    %                               delete(h);
                    %                         end
                    obs_zono = local_to_zono(one_obs);
                    obstacle = obs_zono.Z;
                    
                    %Slice, Do it now! 6 , 7 can be sliced later, taking values vx_des_soft;y_des_soft
                    %Need to get rid of the extra zonos that have negative velocity
                    
                    c = Z(obs_dim, 1);
                    G = Z(:, 2:end);
                    
                    for k_idx = 1:length(k_dim)
                        [~, k_col(k_idx)] = find(G(k_dim(k_idx), :) ~= 0); % find "k-sliceable" generators
                    end
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
            k1c = AH.v_array(vd_idx);k1g = AH.zono_full.kvdg;
            k2c = AH.y_array(y_des_idx);k2g = AH.zono_full.kyg;
            c_k = [k1c; k2c];
            g_k = [k1g; k2g];
            k1_user= vx_des_soft;k2_user=y_des_soft;
            
            % grid over parameter space
            if AH.flags.discrete_flag == false
                k1_sample = linspace(k1c - k1g, k1c + k1g, 13);
                k2_sample = linspace(k2c - k2g, k2c + k2g, 13);
                k1_usr_idx = 5; k2_usr_idx = 5; %start with the center of parameter space
                replace_start_idx = 0;
            else
                k1_sample = linspace(k1c - k1g, k1c + k1g, 3);
                k2_sample = linspace(k2c - k2g, k2c + k2g, 3);
                k1_usr_idx = 2; k2_usr_idx = 2;
                replace_start_idx = 1;
            end
            k1_delta = k1_sample(2)-k1_sample(1);
            k2_delta = k2_sample(2)-k2_sample(1);
            
            avaliable_action_set = zeros(0,3);
            if AH.plot_flag
                figure(9);subplot(1,2,2);
                h1 = scatter([],[],50,[204, 255, 204]/255,'MarkerEdgeColor',[0  0  0],'LineWidth',2)
                % scatter(k1_user, k2_user,100,'y','MarkerEdgeColor',[0  0  0],'LineWidth',2);
            end
            for replace_idx = replace_start_idx:length(k1_sample)*length(k2_sample)
                % was using a sprial search to make algorithm faster.
                % However, the first safe index that the sprial search
                % finds is not the closest one. So instead, just sprial
                % search the whole parameter space and record the avaliable
                % ones with distance
                if replace_idx == 0
                    K = [k1_user; k2_user];
                else
                    [d1,d2] = sprial_index(replace_idx); % search using sprial, deprecated.
                    k1idx = k1_usr_idx + d1;
                    k2idx = k2_usr_idx + d2;
                    if k1idx > length(k1_sample)  || k1idx < 1 || k2idx < 1 ||  k2idx > length(k2_sample)
                        continue;
                    end
                    K = [k1_sample(k1idx); k2_sample(k2idx)];
                    %add noise here to make it learn more robustly
                    %make sure it is centered at K and 3 stdev away is the
                    %next data point. (one idea is to try to make the final
                    %coverage uniform.
                    %also this should be turned off during evaluation.
                    % also don't add noise when doing discrete
                    if ~AH.S.eval && AH.flags.discrete_flag == false
                        sigma = k1_delta/5;
                        K = sigma*randn(2,1) + K;
                    end
                    
                end
                
                %final parameter check
                if K(1) > AH.spdub || K(1) < AH.spdlb || K(2) + agent_state(2) > AH.upper_bd || K(2) + agent_state(2) < AH.lower_bd
                    continue;
                end
                
                
                lambdas = (K - c_k)./g_k; % given a parameter, get coefficients on k_slc_G generators
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
                    if AH.plot_flag
                        h1.XData= [h1.XData K(1)];
                        h1.YData= [h1.YData K(2)];
                    end
                    %save avaliable parameters in a set
                    avaliable_action_set = [avaliable_action_set; K(1) K(2) (k1_user-K(1))^2+(k2_user-K(2))^2];
                    if replace_idx == 0
                        break;% just exit if first selected action already safe
                    end
                end
                
            end
            
            if isempty(avaliable_action_set)
                K = [];
                replaced = 1;
            else
                % find the closest action based on 2nd norm
                [~,best_replace]=min(avaliable_action_set(:,3));
                K = avaliable_action_set(best_replace,1:2)';
                if abs(k1_user-K(1)) < (k1_delta - AH.eps) && abs(k2_user-K(2)) < (k2_delta - AH.eps)
                    replaced = 0;
                else
                    replaced = 1;
                end
                if AH.plot_flag
                    figure(9);subplot(1,2,2)
                    %                     scatter(K(1),K(2),'g');
                    scatter(K(1),K(2),200,'g','filled','MarkerEdgeColor','k','LineWidth',2)
                    
                end
            end
            
        end
        function collision_array = check_zono_collide_halfspace(AH,O,agent_state,v_ini_idx,y_des_idx,h_ini_idx,del_idx,vd_idx)
            % same as find_replace_action execpt check the full parameter range in a bin, for
            % plotting purposes only
            %             if AH.plot_flag
            % %                 figure(2);clf;hold on;axis equal;
            %             end
            zono_peak_peak = AH.zono_full.res{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx,1};
            zono_peak_mid = AH.zono_full.res{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx,2};
            zono_peak_stop = AH.zono_full.res{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx,3};
            error_data_cur = AH.zono_full.error_table{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx};
            %             zono_one_all = [zono_peak_peak; zono_peak_mid;zono_peak_stop];
            
            zono_peak_one = [zono_peak_peak; zono_peak_mid;zono_peak_stop];
            %             zono_stop_one = P.zono_stop.res_stop{v_ini_idx ,y_des_idx,h_ini_idx,del_idx,vd_idx};
            n = length(zono_peak_one); %+length(zono_stop_one);
            
            collision_array = zeros(n, 1);
            if n == 0
                return
            end
            

            dummy_state = agent_state;
            dummy_state(3) = 0; % since the zonotope here doesn't need to be rotated.

            O_pts = world_to_local(dummy_state,O);

            O_pts = O_pts';
            
            
            obs_dim = [1; 2]; % note that the obstacle exists in the x-y space (not theta or v)
            k_dim = [3; 6]; % note that the parameters k are in the 5th and 6th rows of the zonotopes
            buffer_dist = 0; % assume no buffer.
            %
            A_con = {};
            b_con = {};
            s_con = {};
            for t_idx = 1: n
                %                     if t_idx <= length(zono_peak_one)
                zono_one = zono_peak_one{t_idx}{1};

                if isempty(zono_one)
                    error("no such zonotope"+num2str(v_ini_idx)+num2str(y_des_idx)+num2str(h_ini_idx)+num2str(del_idx)+num2str(vd_idx));
                end
                y_des_soft = bound_values(AH.y_des - agent_state(2),AH.zono_full.y_range(y_des_idx)-AH.zono_full.kyg+AH.eps, AH.zono_full.y_range(y_des_idx)+AH.zono_full.kyg-AH.eps);% can't slice on edge
                vx_des_soft =bound_values(AH.vx_des,AH.zono_full.v_range(vd_idx)-AH.zono_full.kvdg+AH.eps, AH.zono_full.v_range(vd_idx)+AH.zono_full.kvdg-AH.eps);
                vi_soft =bound_values(agent_state(4),AH.zono_full.v_range(v_ini_idx)-AH.zono_full.kvig+AH.eps, AH.zono_full.v_range(v_ini_idx)+AH.zono_full.kvig-AH.eps);
                heading_soft = bound_values(agent_state(3), AH.zono_full.h_range(h_ini_idx)-AH.zono_full.h_ini_range+AH.eps, AH.zono_full.h_range(h_ini_idx)+AH.zono_full.h_ini_range-AH.eps);
                zono_one = zonotope_slice(zono_one, [4;5], [vi_soft;heading_soft]);
                
                center_bb = error_data_cur(1:2,t_idx);
                h =  error_data_cur(3,t_idx);
                gen = error_data_cur(4:5,t_idx);
                len = gen(1);
                width = gen(2);
                ego_gen = [[cos(h)*len; sin(h)*len], [sin(-h)*width; cos(-h)*width]];
                gen = zeros(7,2);gen(1:2,1:2) = ego_gen;
                err_zono  = zonotope([[center_bb(1);center_bb(2);0;0;0;0;0], gen]);

                zono_one = (err_zono+zono_one);

                Z = zono_one.Z;
                A_obs_array= []; b_obs_array = [];size_array=[]; size_idx = 0;
                %                  P.vx_des,P.y_des - agent_state(2),agent_state(3),agent_state(4),agent_state(5)
                for obs_idx = 1:(size(O_pts,1)+1)/6
                    one_obs = O_pts((obs_idx-1)*6+1:obs_idx*6-1,:);
                    
                    if max(one_obs(:,1)) < -20 || min(one_obs(:,1)) > 50
                        %                     if(max(one_obs(:,1)) - min(one_obs(:,1)) ) < 100 || max(one_obs(:,2))>0
                        
                        continue;
                    end

                    obs_zono = local_to_zono(one_obs);
                    obstacle = obs_zono.Z;

                    
                    c = Z(obs_dim, 1);
                    G = Z(:, 2:end);
                    
                    for k_idx = 1:length(k_dim)
                        [~, k_col(k_idx)] = find(G(k_dim(k_idx), :) ~= 0); % find "k-sliceable" generators
                    end
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
            k1c = AH.v_array(vd_idx);k1g = AH.zono_full.kvdg;
            k2c = AH.y_array(y_des_idx);k2g = AH.zono_full.kyg;
            c_k = [k1c; k2c];
            g_k = [k1g; k2g];
            k1_user= vx_des_soft;k2_user=y_des_soft;
            lambdas = ([k1_user; k2_user] - c_k)./g_k;
    
            lims = [k1c - k1g, k1c - k1g, k1c + k1g, k1c + k1g, k1c - k1g; k2c - k2g, k2c + k2g, k2c + k2g, k2c - k2g, k2c - k2g];
            
            % grid over parameter space
            k1_sample = linspace(k1c - k1g, k1c + k1g, 20);
            k2_sample = linspace(k2c - k2g, k2c + k2g, 20);
            [Xk, Yk] = meshgrid(k1_sample, k2_sample);
            Zk = inf*ones(size(Xk));
            for i = 1:length(k1_sample)
                for j = 1:length(k2_sample)
                    K = [Xk(i, j); Yk(i, j)];
                    lambdas = (K - c_k)./g_k; % given a parameter, get coefficients on k_slc_G generators
                    for k = 1:length(A_con)
                        Zk_tmp = A_con{k}*lambdas - b_con{k}; % A*lambda - b <= 0 means inside unsafe set
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
            if AH.plot_flag
                % plot parameters and unsafe set
                figure(9); clf;
                set(gcf, 'Position', [1 1 1000 300]);
                subplot(1,2,2);hold on;
                %                 title('', 'FontSize', 15);
                xlabel('$\kappa_{\mathrm{pk}} \mathrm{[m/s]}$', 'Interpreter', 'latex', 'FontSize', 15);
                ylabel('$\kappa_{\mathrm{y}} \mathrm{[m]}$', 'Interpreter', 'latex', 'FontSize', 15);
                
                contourf(Xk, Yk, -Zk, [0, 0], 'FaceColor', [255 51 51]/255); % show zero level set contours
                
                %% have user select a point, see the corresponding slice of FRS:
                
                k1_user= vx_des_soft;k2_user=y_des_soft;
                
                scatter(k1_user, k2_user,200,'y','filled','MarkerEdgeColor',[0  0  0],'LineWidth',2);
                figure(1)
    
            end
        end
    end
end



