Z_des = A.Z_des ;

nan_log = isnan(Z_des(1,:)) ;

N_trajs = sum(nan_log) ;

traj_start_idxs = 1:size(Z_des,2) ;
traj_start_idxs = traj_start_idxs(nan_log) ;

v_peak = [] ;

for idx = 1:(length(traj_start_idxs)-1)
    if idx == 1
        start_idx = 1 ;
        end_idx = traj_start_idxs(1) ;
    else
        start_idx = traj_start_idxs(idx) + 1 ;
        end_idx = traj_start_idxs(idx+1) - 1 ;
    end
    
    Z_des_idx = Z_des(:,start_idx:end_idx) ;
    
    v_peak = [v_peak, max(Z_des_idx(4:6,:),[],2)] ;
end