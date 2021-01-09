function [E_lo, E_hi] = put_error_in_time_bins_2d(E,T,t_lo,t_hi)
    % E is a 3-by-Nt of error samples relative to the trajectories
    % T is a 1-by-NT of times corresponding to the error samples
    % t_lo is an N-by-1 of the lower times of each error time bin
    % t_hi is an N-by-1 of the higher times of each error time bin
    %
    % E_lo is an N-by-3 of the (x,y,z) error lower bounds in each time bin
    % defined by t_lo and t_hi, and E_hi is an N-by-3 of the higher error
    % bounds similarly
    %
    % Testing putting error in time bins:
    % t_lo = 0:0.01:(3-0.01) ;
    % t_hi = 0.01:0.01:3 ;
    % T = 0:0.001:3 ;
    % E = 2*rand(3,length(T)) - 1;
    % [E_lo, E_hi] = put_error_in_time_bins(E,T,t_lo(:),t_hi(:)) ;
    
    % create matrix to logically index time bins
    Nt = length(T) ;
    N = length(t_lo) ;
    
    T_mat = repmat(T,N,1) ;
    t_lo_mat = repmat(t_lo,1,Nt) ;
    t_hi_mat = repmat(t_hi,1,Nt) ;
    
    T_bins_log = (T_mat >= t_lo_mat) & (T_mat <= t_hi_mat) ;
    
    % create matrices of x,y,z error
    Ex = repmat(E(1,:),N,1) ;
    Ey = repmat(E(2,:),N,1) ;
%     Ez = repmat(E(3,:),N,1) ;
    
    % nan out error values not in the time bins
    Ex(~T_bins_log) = nan ;
    Ey(~T_bins_log) = nan ;
%     Ez(~T_bins_log) = nan ;
    
    % compute E_lo
    E_lo = [min(Ex,[],2,'omitnan'), min(Ey,[],2,'omitnan')] ;
    E_hi = [max(Ex,[],2,'omitnan'), max(Ey,[],2,'omitnan')] ;
end