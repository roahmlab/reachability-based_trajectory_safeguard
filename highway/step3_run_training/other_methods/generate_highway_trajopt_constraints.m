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
        %split each obstacle
        
        if max(one_obs(:,1)) < -20 || min(one_obs(:,1)) > 50
            % donn't consider obstacles that are too far away
            
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
    s_con{t_idx} = size_array;% since size may not be same for each obstacle constraint, A*lamda can be multiplied togeher, but need to keep result in a size vector, to determine each obstacle intersection seperately.
end

end
