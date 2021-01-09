function plot_filled_dynamic_obstacles(O_cur,start_color,end_color,edge_color,fig_num,buffer,T,face_alpha,direction)
    % parse inputs
    if nargin<9
        direction = 'ascend';
    if nargin <8
        face_alpha = 1.0;
    if nargin <7
        T=Inf;
    if nargin <6
        buffer = 0;
    if nargin < 5
        fig_num = 1 ;
    if nargin < 4
        edge_color = 'none' ;
    if nargin < 3
        end_color = [1 0.3 0.3] ;
    if nargin < 2
        start_color = [1 0.9 0.9] ;
    end
    end
    end
    end
    end
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
    t_plot = sort(t_plot,direction) ;
    t_max = max(t_plot) ;
    t_min = min(t_plot) ;
    
    % for each time, plot patches of the relevant obstacles
    for t_idx = t_plot
        if t_idx<=T
        t_log = O_cur(3,:) == t_idx; %| isnan(O_cur(1,:)) ;
        O_idx = O_cur(:,t_log) ;
        
        N_polys = sum(isnan(O_idx(1,:)));
        f_polys = find(isnan(O_idx(1,:)));
        f_polys = [0,f_polys];
%         X_data = reshape(O_idx(1,:),[],N_polys);
%                
%         Y_data = reshape(O_idx(2,:),[],N_polys);
%         
%         %filter out nans
%         X_data = reshape(X_data(~isnan(X_data)),[],N_polys);
%         
%         Y_data = reshape(Y_data(~isnan(Y_data)),[],N_polys);
        
        if t_max==t_min || all(end_color == start_color)
            color_cur = start_color;
        else
        
         color_cur = start_color + (t_idx / (t_max - t_min)).*(end_color - start_color) ;
        end
        

        for i=1:N_polys
            Pbuff = bufferPolygonObstacles(O_idx(1:2,f_polys(i)+1:f_polys(i+1)-1),buffer,'square');
            X_buff=Pbuff(1,:)';
            Y_buff=Pbuff(2,:)';
            patch(X_buff,Y_buff,color_cur,'EdgeColor',edge_color,'FaceAlpha',face_alpha)
        end
      


        end
    end
end