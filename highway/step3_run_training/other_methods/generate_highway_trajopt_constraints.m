function [A_con, b_con, s_con] = generate_highway_trajopt_constraints(agent_state,...
    FRS, O_pts)

obs_dim = [1; 2];
k_dim =[3; 6];
%zono_one = zonotope_slice(zono_one, [4;5], [vi_soft;h_soft]);
IC_dim = [4; 5];
IC = [agent_state(4);agent_state(3);] ;
buffer_dist = 0;
A_con = {};
b_con = {};
s_con = {};

%% compute the unsafe set for each FRS time slice
for t_idx = 1:length(FRS)
    frs = FRS{t_idx};
    
    % slice frs by initial condition
    frs = zonotope_slice(frs, IC_dim, IC);
    
    Z = frs.Z;
    A_obs_array= []; b_obs_array = [];size_array=[]; size_idx = 0;
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

end

% function [R_unsafe] = compute_unsafe_params(frs, obs, workspace_dim, param_dim)
% % frs is a zonotope representing the FRS at the current time
% % obs is a cell containing obstacles represented as zonotopes
% % workspace_dim is a vector containing the indices of the workspace
% % dimensions in the FRS
% % param_dim is a vector containing the indices of the workspace dimensions
% % in the FRS.
%
% frs_Z = get(frs, 'Z');
% frs_c = frs_Z(:, 1);
% frs_G = frs_Z(:, 2:end);
%
% % get the indices of generators corresponding to parameters
% frs_param_idx = [];
% for idx = 1:length(param_dim)
%     myidx = find(frs_G(param_dim(idx), :) ~= 0);
%     switch length(myidx)
%         case 0
%             error('No generator for param_dim %i', param_dim(idx));
%         case 1
%             % all good
%         otherwise
%             error('More than one generator for param_dim %i', param_dim(idx));
%     end
%     frs_param_idx(idx, 1) = myidx;
% end
%
% % separate out the "unparameterized" part of FRS
% frs_unparameterized_c = zeros(length(workspace_dim), 1);
% frs_unparameterized_G = frs_G(workspace_dim, :);
% frs_unparameterized_G(:, frs_param_idx) = [];
%
% % get unsafe params for each obstacle
% R_unsafe = {};
% for idx = 1:length(obs)
%     % buffer obstacle by unparameterized portion of FRS
%     % this is ok, since our dynamics are independent of position!
%     buffered_obs_c = frs_unparameterized_c + obs{idx}.Z(:, 1);
%     buffered_obs_G = [frs_unparameterized_G, obs{idx}.Z(:, 2:end)];
%     % subtract off current center position
%     buffered_obs_c = buffered_obs_c - frs_c(workspace_dim, 1);
%     % turn buffered obstacle zonotope into an interval.
%
%     % determine left and right limit
%     delta = sum(abs(buffered_obs_G),2);
%     leftLimit = buffered_obs_c - delta;
%     rightLimit = buffered_obs_c + delta;
%
%     unsafe_c = zeros(length(param_dim), 1);
%     unsafe_G = zeros(length(param_dim), length(param_dim));
%
%     for j = 1:length(param_dim)
%         myg = frs_G(:, frs_param_idx(j));
%         mydim = myg(workspace_dim, 1) ~= 0; % get workspace dim this param is associate with.
%
%         beta_1 = leftLimit(mydim)/myg(workspace_dim(mydim), 1);
%         beta_2 = rightLimit(mydim)/myg(workspace_dim(mydim), 1);
%
%         if beta_1 > 1 || beta_2 < -1
%             % no intersection possible.
%             beta_1 = nan;
%             beta_2 = nan;
%         else
%             beta_1 = min(1, max(-1, beta_1));
%             beta_2 = min(1, max(-1, beta_2));
%         end
%         myg_param = frs_G(param_dim, frs_param_idx(j));
%         unsafe_c = unsafe_c + (beta_1 + beta_2)/2*myg_param;
%         unsafe_G(:, j) = (beta_2 - beta_1)/2*myg_param;
%     end
%     % add back on original center of zono
%     unsafe_c = unsafe_c + frs_c(param_dim);
%     unsafe_Z = [unsafe_c, unsafe_G];
%     if any(any(isnan(unsafe_Z)))
%         R_unsafe{idx} = {};
%     else
%         for k = 1:length(param_dim)
%             plus_g = unsafe_c(k) + unsafe_G(k, k);
%             minus_g = unsafe_c(k) - unsafe_G(k, k);
%             R_unsafe{idx}(k, 1) = min(plus_g, minus_g);
%             R_unsafe{idx}(k, 2) = max(plus_g, minus_g);
%         end
%     end
% end
% end