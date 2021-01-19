function [A_con, b_con, s_con] = generate_cartpole_trajopt_constraints(AH,agent_state,...
    ka_idx,v_idx,kpkidx,theta_idx,theta_dot_idx, O)
    
    vd_idx_relative = kpkidx-v_idx+2;
    zono_peak     = AH.zono_full.res{v_idx,ka_idx,vd_idx_relative,1};
    zono_failsafe = AH.zono_full.res{v_idx,ka_idx,vd_idx_relative,2};
    error_data_cur = AH.zono_full.error_table{v_idx,ka_idx,vd_idx_relative,theta_idx,theta_dot_idx};
    zono_prs = [zono_peak; zono_failsafe];
    n = length(zono_prs);
    if n == 0
        warning("no such zonotope kv,ka,kpeak ="+num2str(v_idx)+num2str(ka_idx)+num2str(vd_idx));
        return
    end
%agent state already bounded in above step
%     kv_soft = bound_values(agent_state(2),AH.zono_full.vi(v_idx)-AH.zono_full.kvg,AH.zono_full.vi(v_idx)+AH.zono_full.kvg);
%     ka_soft = bound_values(agent_state(5),AH.zono_full.ka(ka_idx)-AH.zono_full.kag,AH.zono_full.ka(ka_idx)+AH.zono_full.kag);
    IC = [agent_state(2);agent_state(5)];
    IC_dim = [3;4];

    pos_cart = agent_state(1);
    O_pts = O - pos_cart;
    
    obst_dim = [1];%integrated velocity = x
    k_dim = [2];% peak or desired velocity dimension
    A_con = {};
    b_con = {};
    s_con = {};
    buffer_dist = 0;

    %% compute the unsafe set for each FRS time slice
    for t_idx = 1: n
            zono_one = zono_prs{t_idx}{1};
%             kv_soft = bound_values(agent_state(2),AH.zono_full.vi(v_idx)-AH.zono_full.kvg,AH.zono_full.vi(v_idx)+AH.zono_full.kvg);
%             ka_soft = bound_values(agent_state(5),AH.zono_full.ka(ka_idx)-AH.zono_full.kag,AH.zono_full.ka(ka_idx)+AH.zono_full.kag);
%                 %                 kvd_soft =    
            zono_one = zonotope_slice(zono_one, IC_dim, IC);
            centererr = mean(error_data_cur(t_idx,:));
            gen_err = zeros(5,5); gen_err(1,1) = (error_data_cur(t_idx,2)- error_data_cur(t_idx,1))/2;
            err_zono  = zonotope([[centererr;0;0;0;0], gen_err]);
                
            zono_one = (err_zono+zono_one);

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
end
