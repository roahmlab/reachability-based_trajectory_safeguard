function plot_gradient_dynamic_points(O_cur,start_color,end_color,fig_num)
    % parse inputs
    if nargin < 4
        fig_num = 1 ;
    if nargin < 3
        end_color = [1 0.3 0.3] ;
    if nargin < 2
        start_color = [1 0.9 0.9] ;
    end
    end
    end

    
    % set up figure
    figure(fig_num) ; hold on ;
    
    % clean up nans in O_cur; add a column of nans at the end to help
    % indexing
     if isnan(O_cur(1,1))
        O_cur = O_cur(:,2:end) ;
    end
    
    if ~isnan(O_cur(1,end))
        O_cur = [O_cur, nan(3,1)] ;
    end
    
    % get times to plot
    t_plot = unique(O_cur(3,:)) ;
    t_plot = t_plot(~isnan(t_plot)) ;
    t_plot = sort(t_plot,'ascend') ;
    t_max = max(t_plot) ;
    t_min = min(t_plot) ;
    
    % for each time, plot patches of the relevant obstacles
    for t_idx = t_plot
        t_log = O_cur(3,:) == t_idx; %| isnan(O_cur(1,:)) ;
        O_idx = O_cur(:,t_log) ;
        

        
        color_cur = start_color + (t_idx / (t_max - t_min)).*(end_color - start_color) ;
        
         plot(O_idx(1,:),O_idx(2,:),'.','Color',color_cur)

    end
end