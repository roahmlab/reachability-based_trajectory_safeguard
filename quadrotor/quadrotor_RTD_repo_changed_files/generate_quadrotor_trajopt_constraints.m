function [A_con, b_con] = generate_quadrotor_trajopt_constraints(v_0, a_0,...
    FRS, obstacles, tracking_error_in)

    position_dimensions = [1; 6; 11];
    param_dimensions = [4; 9; 14];
    IC_dim = [2; 7; 12; 3; 8; 13];
    IC = [v_0; a_0] ;
    
    %% vars used often
    nargin_check = nargin > 4 ;

    if nargin_check
        tracking_error_in_class = class(tracking_error_in) ;
        switch tracking_error_in_class
            case 'table'
                tbl = get_tracking_error_for_current_velocity(tracking_error_in, v_0);
            case 'zonotope'
                zono = tracking_error_in ;
        end
    end

    %% compute the unsafe set for each FRS time slice
    b_con = [];
    for idx = 1:length(FRS)
        frs = FRS{idx}{1};

        % slice frs by initial condition
%         try
            frs = zonotope_slice(frs, IC_dim, IC);
%         catch
%             dbstop in generate_quadrotor_trajopt_constraints at 32
%              disp('hi2')
%         end
        
        
        % add in tracking error
        if nargin_check
            switch tracking_error_in_class
                case 'table'
                    try
                        frs_with_error = frs + tbl(idx, :).zonotope{1};
                    catch
                        dbstop in generate_quadrotor_trajopt_constraints at 44
                        disp('hi')
                    end
                case 'zonotope'
                    frs_with_error = frs + zono ;
            end
        else
            frs_with_error = frs ;
        end

        % compute unsafe trajectory parameters
        K_unsafe = compute_unsafe_params(frs_with_error, obstacles{idx},...
                        position_dimensions, param_dimensions);
        
        % create the constraints
        % unsafeK is a 3x2 matrix, with the lower bound and upper bound of
        % unsafeK in each axis.
        for j = 1:length(K_unsafe)
            if ~isempty(K_unsafe{j})
                b = [-K_unsafe{j}(1, 1); K_unsafe{j}(1, 2);
                    -K_unsafe{j}(2, 1); K_unsafe{j}(2, 2);
                    -K_unsafe{j}(3, 1); K_unsafe{j}(3, 2)];
                b_con = [b_con; b];
            end
        end
    end

    A = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1] ;
    A_con = repmat(A,size(b_con,1)/6,1) ;
end

function [R_unsafe] = compute_unsafe_params(frs, obs, workspace_dim, param_dim)
% frs is a zonotope representing the FRS at the current time
% obs is a cell containing obstacles represented as zonotopes
% workspace_dim is a vector containing the indices of the workspace
% dimensions in the FRS
% param_dim is a vector containing the indices of the workspace dimensions
% in the FRS.

    frs_Z = get(frs, 'Z');
    frs_c = frs_Z(:, 1);
    frs_G = frs_Z(:, 2:end);

    % get the indices of generators corresponding to parameters
    frs_param_idx = [];
    for idx = 1:length(param_dim)
        myidx = find(frs_G(param_dim(idx), :) ~= 0);
        switch length(myidx)
            case 0
                error('No generator for param_dim %i', param_dim(idx));
            case 1
                % all good
            otherwise
                error('More than one generator for param_dim %i', param_dim(idx));
        end
        frs_param_idx(idx, 1) = myidx;
    end

    % separate out the "unparameterized" part of FRS
    frs_unparameterized_c = zeros(length(workspace_dim), 1);
    frs_unparameterized_G = frs_G(workspace_dim, :);
    frs_unparameterized_G(:, frs_param_idx) = [];

    % get unsafe params for each obstacle
    R_unsafe = {};
    for idx = 1:length(obs)
        % buffer obstacle by unparameterized portion of FRS
        % this is ok, since our dynamics are independent of position!
        buffered_obs_c = frs_unparameterized_c + obs{idx}.Z(:, 1);
        buffered_obs_G = [frs_unparameterized_G, obs{idx}.Z(:, 2:end)];
        % subtract off current center position
        buffered_obs_c = buffered_obs_c - frs_c(workspace_dim, 1);
        % turn buffered obstacle zonotope into an interval.
        
        % determine left and right limit
        delta = sum(abs(buffered_obs_G),2);
        leftLimit = buffered_obs_c - delta;
        rightLimit = buffered_obs_c + delta;

        unsafe_c = zeros(length(param_dim), 1);
        unsafe_G = zeros(length(param_dim), length(param_dim));
        
        for j = 1:length(param_dim)
            myg = frs_G(:, frs_param_idx(j));
            mydim = myg(workspace_dim, 1) ~= 0; % get workspace dim this param is associate with.

            beta_1 = leftLimit(mydim)/myg(workspace_dim(mydim), 1);
            beta_2 = rightLimit(mydim)/myg(workspace_dim(mydim), 1);

            if beta_1 > 1 || beta_2 < -1
                % no intersection possible.
                beta_1 = nan;
                beta_2 = nan;
            else
                beta_1 = min(1, max(-1, beta_1));
                beta_2 = min(1, max(-1, beta_2));
            end
            myg_param = frs_G(param_dim, frs_param_idx(j));
            unsafe_c = unsafe_c + (beta_1 + beta_2)/2*myg_param;
            unsafe_G(:, j) = (beta_2 - beta_1)/2*myg_param;
        end
        % add back on original center of zono
        unsafe_c = unsafe_c + frs_c(param_dim);
        unsafe_Z = [unsafe_c, unsafe_G];
        if any(any(isnan(unsafe_Z)))
            R_unsafe{idx} = {};
        else
            for k = 1:length(param_dim)
                plus_g = unsafe_c(k) + unsafe_G(k, k);
                minus_g = unsafe_c(k) - unsafe_G(k, k);
                R_unsafe{idx}(k, 1) = min(plus_g, minus_g);
                R_unsafe{idx}(k, 2) = max(plus_g, minus_g);
            end
        end
    end
end
