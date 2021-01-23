classdef CartPoleVisualizer < rl.env.viz.AbstractFigureVisualizer
% CARTPOLEVISUALIZER

% Copyright 2018 The MathWorks, Inc.

    methods
        function this = CartPoleVisualizer(A)
            this = this@rl.env.viz.AbstractFigureVisualizer(A);
        end
    end
    methods (Access = protected)
        function f = buildFigure(this)
            
            f = figure(...
                'Toolbar','none',...
                'Visible','on',...
                'HandleVisibility','off', ...
                'NumberTitle','off',...
                'Name',getString(message('rl:env:VizNameCartPole')),... 
                'MenuBar','none',...
                'CloseRequestFcn',@(~,~)delete(this));
            if ~strcmp(f.WindowStyle,'docked')
                f.Position(3:4) = [600 200];
            end
            ha = gca(f);
            
            ha.XLimMode = 'manual';
            ha.YLimMode = 'manual';
            ha.ZLimMode = 'manual';
            ha.DataAspectRatioMode = 'manual';
            ha.PlotBoxAspectRatioMode = 'manual';
            ha.YTick = [];
            
            ha.XLim = [-10 10]; %[-5 5]
            ha.YLim = [-3 3];    %[0 3];
            
            hold(ha,'on');
        end
        function updatePlot(this)

            %env = this.Environment;
            A = this.A;
                
            f = this.Figure;
            ha = gca(f);
            
            % ch is the height of the cart
            ch = 0.4;
            % cw is the width of the cart
            cw = 0.4;
            % ph is the length of the pole
            ph = A.Length*2;
            % pw is the width of the pole
            pw = cw/2;
            % theta limits
%             theta_lim = env.ThetaThresholdRadians;
%             y_theta_lim_ = cos(theta_lim)*ph*1.5;
%             x_theta_lim_ = sin(theta_lim)*ph*1.5;
            x_x_lim = 3.5;
            y_x_lim = max(ha.YLim);
            
            % extract position and angle
            state = A.AH.A.state(:,end);
            x = state(1);
            theta = -state(3);%-state(3);
            
            cartplot = findobj(ha,'Tag','cartplot');
            poleplot = findobj(ha,'Tag','poleplot');
            
            leftblockplot = findobj(ha,'Tag','leftblockplot');
            rightblockplot = findobj(ha,'Tag','rightblockplot');
            
            theta_lim_lt_plot = findobj(ha,'Tag','thetalimltplot');
            theta_lim_rt_plot = findobj(ha,'Tag','thetalimrtplot');
            
            
            lower_traj = this.A.AH.lowerbound_array;
            higher_traj = this.A.AH.higherbound_array;
            
            n = length(higher_traj);
            
            traj_enable = 0;
            
            if traj_enable
            safe_estimate_plot = findobj(ha,'Tag','safe_estimate_plot');
            unsafe_estimate_plot = findobj(ha,'Tag','unsafe_estimate_plot');
            safe_estimate_plot_2 = findobj(ha,'Tag','safe_estimate_plot_2');
            unsafe_estimate_plot_2 = findobj(ha,'Tag','unsafe_estimate_plot_2');
            safe_estimate_plot_3 = findobj(ha,'Tag','safe_estimate_plot_3');
            unsafe_estimate_plot_3 = findobj(ha,'Tag','unsafe_estimate_plot_3');
            safe_estimate_plot_4 = findobj(ha,'Tag','safe_estimate_plot_4');
            unsafe_estimate_plot_4 = findobj(ha,'Tag','unsafe_estimate_plot_4');
            safe_estimate_plot_5 = findobj(ha,'Tag','safe_estimate_plot_5');
            unsafe_estimate_plot_5 = findobj(ha,'Tag','unsafe_estimate_plot_5');
            safe_estimate_plot_6 = findobj(ha,'Tag','safe_estimate_plot_6');
            unsafe_estimate_plot_6 = findobj(ha,'Tag','unsafe_estimate_plot_6');
            safe_estimate_plot_7 = findobj(ha,'Tag','safe_estimate_plot_7');
            unsafe_estimate_plot_7 = findobj(ha,'Tag','unsafe_estimate_plot_7');
            safe_estimate_plot_8 = findobj(ha,'Tag','safe_estimate_plot_8');
            unsafe_estimate_plot_8 = findobj(ha,'Tag','unsafe_estimate_plot_8');
            safe_estimate_plot_9 = findobj(ha,'Tag','safe_estimate_plot_9');
            unsafe_estimate_plot_9 = findobj(ha,'Tag','unsafe_estimate_plot_9');
            safe_estimate_plot_10 = findobj(ha,'Tag','safe_estimate_plot_10');
            unsafe_estimate_plot_10 = findobj(ha,'Tag','unsafe_estimate_plot_10');
            unsafe_estimate_plot_11 = findobj(ha,'Tag','unsafe_estimate_plot_11');
            unsafe_estimate_plot_12 = findobj(ha,'Tag','unsafe_estimate_plot_12');
            unsafe_estimate_plot_13 = findobj(ha,'Tag','unsafe_estimate_plot_13');
            unsafe_estimate_plot_14 = findobj(ha,'Tag','unsafe_estimate_plot_14');
            unsafe_estimate_plot_15 = findobj(ha,'Tag','unsafe_estimate_plot_15');
            unsafe_estimate_plot_16 = findobj(ha,'Tag','unsafe_estimate_plot_16');
            unsafe_estimate_plot_17 = findobj(ha,'Tag','unsafe_estimate_plot_17');
            unsafe_estimate_plot_18 = findobj(ha,'Tag','unsafe_estimate_plot_18');
            unsafe_estimate_plot_19 = findobj(ha,'Tag','unsafe_estimate_plot_19');
            unsafe_estimate_plot_20 = findobj(ha,'Tag','unsafe_estimate_plot_20');
            end
            
            x_lim_lt_plot = findobj(ha,'Tag','xlimltplot');
            x_lim_rt_plot = findobj(ha,'Tag','xlimrtplot');
            
            if traj_enable
            
            if isempty(cartplot) || ~isvalid(cartplot) ...
                    || isempty(poleplot) || ~isvalid(poleplot) ...
                    || isempty(theta_lim_lt_plot)       || ~isvalid(theta_lim_lt_plot)       ...
                    || isempty(theta_lim_rt_plot)       || ~isvalid(theta_lim_rt_plot)       ...
                    || isempty(x_lim_lt_plot)           || ~isvalid(x_lim_lt_plot)           ...
                    || isempty(x_lim_rt_plot)           || ~isvalid(x_lim_rt_plot)           ...
                    || isempty(leftblockplot)           || ~isvalid(leftblockplot)           ...
                    || isempty(rightblockplot)          || ~isvalid(rightblockplot)          ...
                    || isempty(safe_estimate_plot)   || ~isvalid(safe_estimate_plot)            ...
                    || isempty(unsafe_estimate_plot) || ~isvalid(unsafe_estimate_plot)          ...
                    || isempty(safe_estimate_plot_2)   || ~isvalid(safe_estimate_plot_2)        ...
                    || isempty(unsafe_estimate_plot_2) || ~isvalid(unsafe_estimate_plot_2)      ...
                    || isempty(safe_estimate_plot_3)   || ~isvalid(safe_estimate_plot_3)        ...
                    || isempty(unsafe_estimate_plot_3) || ~isvalid(unsafe_estimate_plot_3)      ...
                    || isempty(safe_estimate_plot_4)   || ~isvalid(safe_estimate_plot_4)        ...
                    || isempty(unsafe_estimate_plot_4) || ~isvalid(unsafe_estimate_plot_4)      ...
                    || isempty(safe_estimate_plot_5)   || ~isvalid(safe_estimate_plot_5)        ...
                    || isempty(unsafe_estimate_plot_5) || ~isvalid(unsafe_estimate_plot_5)      ...
                    || isempty(safe_estimate_plot_6)   || ~isvalid(safe_estimate_plot_6)        ...
                    || isempty(unsafe_estimate_plot_6) || ~isvalid(unsafe_estimate_plot_6) ...
                    || isempty(safe_estimate_plot_7)   || ~isvalid(safe_estimate_plot_7)        ...
                    || isempty(unsafe_estimate_plot_7) || ~isvalid(unsafe_estimate_plot_7) ...
                    || isempty(safe_estimate_plot_8)   || ~isvalid(safe_estimate_plot_8)        ...
                    || isempty(unsafe_estimate_plot_8) || ~isvalid(unsafe_estimate_plot_8) ...
                    || isempty(safe_estimate_plot_9)   || ~isvalid(safe_estimate_plot_9)        ...
                    || isempty(unsafe_estimate_plot_9) || ~isvalid(unsafe_estimate_plot_9) ...
                    || isempty(safe_estimate_plot_10)   || ~isvalid(safe_estimate_plot_10)      ...
                    || isempty(unsafe_estimate_plot_10) || ~isvalid(unsafe_estimate_plot_10) ...
                    || isempty(unsafe_estimate_plot_11) || ~isvalid(unsafe_estimate_plot_11) ...
                    || isempty(unsafe_estimate_plot_12) || ~isvalid(unsafe_estimate_plot_12) ...
                    || isempty(unsafe_estimate_plot_13) || ~isvalid(unsafe_estimate_plot_13) ...
                    || isempty(unsafe_estimate_plot_14) || ~isvalid(unsafe_estimate_plot_14) ...
                    || isempty(unsafe_estimate_plot_15) || ~isvalid(unsafe_estimate_plot_15) ...
                    || isempty(unsafe_estimate_plot_16) || ~isvalid(unsafe_estimate_plot_16) ...
                    || isempty(unsafe_estimate_plot_17) || ~isvalid(unsafe_estimate_plot_17) ...
                    || isempty(unsafe_estimate_plot_18) || ~isvalid(unsafe_estimate_plot_18) ...
                    || isempty(unsafe_estimate_plot_19) || ~isvalid(unsafe_estimate_plot_19) ...
                    || isempty(unsafe_estimate_plot_20) || ~isvalid(unsafe_estimate_plot_20) ...

                delete(cartplot);
                delete(poleplot);
                delete(theta_lim_lt_plot);
                delete(theta_lim_rt_plot);
                delete(x_lim_lt_plot);
                delete(x_lim_rt_plot);
                delete(leftblockplot);
                delete(rightblockplot);
                
                delete(safe_estimate_plot);
                delete(safe_estimate_plot_2);
                delete(safe_estimate_plot_3);
                delete(safe_estimate_plot_4);
                delete(safe_estimate_plot_5);
                delete(safe_estimate_plot_6);
                delete(safe_estimate_plot_7);
                delete(safe_estimate_plot_8);
                delete(safe_estimate_plot_9);
                delete(safe_estimate_plot_10);

                delete(unsafe_estimate_plot);
                delete(unsafe_estimate_plot_2);
                delete(unsafe_estimate_plot_3);
                delete(unsafe_estimate_plot_4);
                delete(unsafe_estimate_plot_5);
                delete(unsafe_estimate_plot_6);
                delete(unsafe_estimate_plot_7);
                delete(unsafe_estimate_plot_8);
                delete(unsafe_estimate_plot_9);
                delete(unsafe_estimate_plot_10);
                delete(unsafe_estimate_plot_11);
                delete(unsafe_estimate_plot_12);
                delete(unsafe_estimate_plot_13);
                delete(unsafe_estimate_plot_14);
                delete(unsafe_estimate_plot_15);
                delete(unsafe_estimate_plot_16);
                delete(unsafe_estimate_plot_17);
                delete(unsafe_estimate_plot_18);
                delete(unsafe_estimate_plot_19);
                delete(unsafe_estimate_plot_20);
                
            
                

                % Create polyshape objects
                cartpoly = polyshape([-cw/2,-cw/2,cw/2,cw/2],[0   ,ch       ,ch     ,0   ]);
                polepoly = polyshape([-pw/2,-pw/2,pw/2,pw/2],[ch/2,ph + ch/2,ph+ch/2,ch/2]);
                
                
                x_lim_left_block = [-x_x_lim-0.5, -x_x_lim-0.5,  -x_x_lim, -x_x_lim];
                y_lim_left_block = [-y_x_lim,      y_x_lim,      y_x_lim,  -y_x_lim ];
                x_lim_right_block = [x_x_lim+0.5,  x_x_lim+0.5,   x_x_lim,  x_x_lim];
                y_lim_right_block = [-y_x_lim,     y_x_lim,       y_x_lim, -y_x_lim];
            
                if n~=0
                   
                    lower_width = higher_traj(1);
                    higher_width = higher_traj(2);
                    width_array = [lower_width-cw/2-x,lower_width-cw/2-x,higher_width+cw/2-x,higher_width+cw/2-x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    safe_poly = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(3);
                    higher_width = higher_traj(4);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    safe_poly_2 = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(5);
                    higher_width = higher_traj(6);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    safe_poly_3 = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(7);
                    higher_width = higher_traj(8);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    safe_poly_4 = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(9);
                    higher_width = higher_traj(10);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    safe_poly_5 = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(11);
                    higher_width = higher_traj(12);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    safe_poly_6 = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(13);
                    higher_width = higher_traj(14);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    safe_poly_7 = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(15);
                    higher_width = higher_traj(16);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    safe_poly_8 = polyshape(width_array,height_array);
                                        
                    lower_width = higher_traj(17);
                    higher_width = higher_traj(18);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    safe_poly_9 = polyshape(width_array,height_array);
                                        
                    lower_width = higher_traj(19);
                    higher_width = higher_traj(20);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    safe_poly_10 = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(21);
                    higher_width = higher_traj(22);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(23);
                    higher_width = higher_traj(24);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_2 = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(25);
                    higher_width = higher_traj(26);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_3 = polyshape(width_array,height_array);
                                        
                    lower_width = higher_traj(27);
                    higher_width = higher_traj(28);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_4 = polyshape(width_array,height_array);
                                        
                    lower_width = higher_traj(29);
                    higher_width = higher_traj(30);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_5 = polyshape(width_array,height_array);
                                        
                    lower_width = higher_traj(31);
                    higher_width = higher_traj(32);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_6 = polyshape(width_array,height_array);
                                        
                    lower_width = higher_traj(33);
                    higher_width = higher_traj(34);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_7 = polyshape(width_array,height_array);
                                        
                    lower_width = higher_traj(35);
                    higher_width = higher_traj(36);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_8 = polyshape(width_array,height_array);
                                        
                    lower_width = higher_traj(37);
                    higher_width = higher_traj(38);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_9 = polyshape(width_array,height_array);
                                        
                    lower_width = higher_traj(39);
                    higher_width = higher_traj(40);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_10 = polyshape(width_array,height_array);
                                        
                    lower_width = higher_traj(41);
                    higher_width = higher_traj(42);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_11 = polyshape(width_array,height_array);
                                        
                    lower_width = higher_traj(43);
                    higher_width = higher_traj(44);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_12 = polyshape(width_array,height_array);
                                        
                    lower_width = higher_traj(45);
                    higher_width = higher_traj(46);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_13 = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(47);
                    higher_width = higher_traj(48);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_14 = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(49);
                    higher_width = higher_traj(50);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_15 = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(51);
                    higher_width = higher_traj(52);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_16 = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(53);
                    higher_width = higher_traj(54);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_17 = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(55);
                    higher_width = higher_traj(56);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_18 = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(57);
                    higher_width = higher_traj(58);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_19 = polyshape(width_array,height_array);
                    
                    lower_width = higher_traj(59);
                    higher_width = higher_traj(60);
                    width_array = [lower_width-cw/2 - x,lower_width-cw/2 - x,higher_width+cw/2 - x,higher_width+cw/2 - x];
                    height_array = [-ch,      ch,      ch,  -ch ];
                    unsafe_poly_20 = polyshape(width_array,height_array);
                       
                else
                end

                
                leftblockpoly = polyshape(x_lim_left_block,y_lim_left_block);
                rightblockpoly = polyshape(x_lim_right_block,y_lim_right_block);
                
                % Create the theta limit lines
                lineargs = {...
                    'LineStyle'     ,'--', ...
                    'Color'         ,'k' , ...
                    'LineWidth'     ,1     ...
                    };
%                 theta_lim_lt_plot = line(ha,[0 0],[0 0],lineargs{:});
%                 theta_lim_rt_plot = line(ha,[0 0],[0 0],lineargs{:});
                
                % Create the x limit lines
                x_lim_lt_plot = line(ha,[0 0],[0 0],lineargs{:});
                x_lim_rt_plot = line(ha,[0 0],[0 0],lineargs{:});
                

                cartplot = plot(ha,cartpoly,'FaceColor','b');
                cartplot.Tag = 'cartplot';
                poleplot = plot(ha,polepoly,'FaceColor','b');
                poleplot.Tag = 'poleplot';
                
                leftblockplot = plot(ha,leftblockpoly,'FaceColor','r','FaceAlpha',0.5);
                leftblockplot.Tag = 'leftblockplot';
                rightblockplot = plot(ha,rightblockpoly,'FaceColor','r','FaceAlpha',0.5);
                rightblockplot.Tag = 'rightblockplot';
                
                if n~= 0
                    safe_estimate_plot = plot(ha,safe_poly,'FaceColor','g');
                    safe_estimate_plot_2 = plot(ha,safe_poly_2,'FaceColor','g');
                    safe_estimate_plot_3 = plot(ha,safe_poly_3,'FaceColor','g');
                    safe_estimate_plot_4 = plot(ha,safe_poly_4,'FaceColor','g');
                    safe_estimate_plot_5 = plot(ha,safe_poly_5,'FaceColor','g');
                    safe_estimate_plot_6 = plot(ha,safe_poly_6,'FaceColor','g');
                    safe_estimate_plot_7 = plot(ha,safe_poly_7,'FaceColor','g');
                    safe_estimate_plot_8 = plot(ha,safe_poly_8,'FaceColor','g');
                    safe_estimate_plot_9 = plot(ha,safe_poly_9,'FaceColor','g');
                    safe_estimate_plot_10 = plot(ha,safe_poly_10,'FaceColor','g');


                    safe_estimate_plot.Tag = 'safe_estimate_plot';
                    safe_estimate_plot_2.Tag = 'safe_estimate_plot_2';
                    safe_estimate_plot_3.Tag = 'safe_estimate_plot_3';
                    safe_estimate_plot_4.Tag = 'safe_estimate_plot_4';
                    safe_estimate_plot_5.Tag = 'safe_estimate_plot_5';
                    safe_estimate_plot_6.Tag = 'safe_estimate_plot_6';
                    safe_estimate_plot_7.Tag = 'safe_estimate_plot_7';
                    safe_estimate_plot_8.Tag = 'safe_estimate_plot_8';
                    safe_estimate_plot_9.Tag = 'safe_estimate_plot_9';
                    safe_estimate_plot_10.Tag = 'safe_estimate_plot_10';

                    unsafe_estimate_plot = plot(ha,unsafe_poly,'FaceColor','y');
                    unsafe_estimate_plot_2 = plot(ha,unsafe_poly_2,'FaceColor','y');
                    unsafe_estimate_plot_3 = plot(ha,unsafe_poly_3,'FaceColor','y');
                    unsafe_estimate_plot_4 = plot(ha,unsafe_poly_4,'FaceColor','y');
                    unsafe_estimate_plot_5 = plot(ha,unsafe_poly_5,'FaceColor','y');
                    unsafe_estimate_plot_6 = plot(ha,unsafe_poly_6,'FaceColor','y');
                    unsafe_estimate_plot_7 = plot(ha,unsafe_poly_7,'FaceColor','y');
                    unsafe_estimate_plot_8 = plot(ha,unsafe_poly_8,'FaceColor','y');
                    unsafe_estimate_plot_9 = plot(ha,unsafe_poly_9,'FaceColor','y');
                    unsafe_estimate_plot_10 = plot(ha,unsafe_poly_10,'FaceColor','y');
                    unsafe_estimate_plot_11 = plot(ha,unsafe_poly_11,'FaceColor','y');
                    unsafe_estimate_plot_12 = plot(ha,unsafe_poly_12,'FaceColor','y');
                    unsafe_estimate_plot_13 = plot(ha,unsafe_poly_13,'FaceColor','y');
                    unsafe_estimate_plot_14 = plot(ha,unsafe_poly_14,'FaceColor','y');
                    unsafe_estimate_plot_15 = plot(ha,unsafe_poly_15,'FaceColor','y');
                    unsafe_estimate_plot_16 = plot(ha,unsafe_poly_16,'FaceColor','y');
                    unsafe_estimate_plot_17 = plot(ha,unsafe_poly_17,'FaceColor','y');
                    unsafe_estimate_plot_18 = plot(ha,unsafe_poly_18,'FaceColor','y');
                    unsafe_estimate_plot_19 = plot(ha,unsafe_poly_19,'FaceColor','y');
                    unsafe_estimate_plot_20 = plot(ha,unsafe_poly_20,'FaceColor','y');


                    unsafe_estimate_plot.Tag = 'unsafe_estimate_plot';
                    unsafe_estimate_plot_2.Tag = 'unsafe_estimate_plot_2';
                    unsafe_estimate_plot_3.Tag = 'unsafe_estimate_plot_3';
                    unsafe_estimate_plot_4.Tag = 'unsafe_estimate_plot_4';
                    unsafe_estimate_plot_5.Tag = 'unsafe_estimate_plot_5';
                    unsafe_estimate_plot_6.Tag = 'unsafe_estimate_plot_6';
                    unsafe_estimate_plot_7.Tag = 'unsafe_estimate_plot_7';
                    unsafe_estimate_plot_8.Tag = 'unsafe_estimate_plot_8';
                    unsafe_estimate_plot_9.Tag = 'unsafe_estimate_plot_9';
                    unsafe_estimate_plot_10.Tag = 'unsafe_estimate_plot_10';
                    unsafe_estimate_plot_11.Tag = 'unsafe_estimate_plot_11';
                    unsafe_estimate_plot_12.Tag = 'unsafe_estimate_plot_12';
                    unsafe_estimate_plot_13.Tag = 'unsafe_estimate_plot_13';
                    unsafe_estimate_plot_14.Tag = 'unsafe_estimate_plot_14';
                    unsafe_estimate_plot_15.Tag = 'unsafe_estimate_plot_15';
                    unsafe_estimate_plot_16.Tag = 'unsafe_estimate_plot_16';
                    unsafe_estimate_plot_17.Tag = 'unsafe_estimate_plot_17';
                    unsafe_estimate_plot_18.Tag = 'unsafe_estimate_plot_18';
                    unsafe_estimate_plot_19.Tag = 'unsafe_estimate_plot_19';
                    unsafe_estimate_plot_20.Tag = 'unsafe_estimate_plot_20';

                end
                

                
%                 theta_lim_lt_plot.Tag = 'thetalimltplot';
%                 theta_lim_rt_plot.Tag = 'thetalimrtplot';
                
                x_lim_lt_plot.Tag = 'xlimltplot';
                x_lim_rt_plot.Tag = 'xlimrtplot';
            else
                cartpoly = cartplot.Shape;
                polepoly = poleplot.Shape;
                leftblockpoly = leftblockplot.Shape;
                rightblockpoly = rightblockplot.Shape;
                
%                 safe_poly = safe_estimate_plot.Shape;
%                 unsafe_poly = unsafe_estimate_plot.Shape;
                
                safe_poly = safe_estimate_plot.Shape;
                safe_poly_2 = safe_estimate_plot_2.Shape;
                safe_poly_3 = safe_estimate_plot_3.Shape;
                safe_poly_4 = safe_estimate_plot_4.Shape;
                safe_poly_5 = safe_estimate_plot_5.Shape;
                safe_poly_6 = safe_estimate_plot_6.Shape;
                safe_poly_7 = safe_estimate_plot_7.Shape;
                safe_poly_8 = safe_estimate_plot_8.Shape;
                safe_poly_9 = safe_estimate_plot_9.Shape;
                safe_poly_10 = safe_estimate_plot_10.Shape;

                unsafe_poly = unsafe_estimate_plot.Shape;
                unsafe_poly_2 = unsafe_estimate_plot_2.Shape;
                unsafe_poly_3 = unsafe_estimate_plot_3.Shape;
                unsafe_poly_4 = unsafe_estimate_plot_4.Shape;
                unsafe_poly_5 = unsafe_estimate_plot_5.Shape;
                unsafe_poly_6 = unsafe_estimate_plot_6.Shape;
                unsafe_poly_7 = unsafe_estimate_plot_7.Shape;
                unsafe_poly_8 = unsafe_estimate_plot_8.Shape;
                unsafe_poly_9 = unsafe_estimate_plot_9.Shape;
                unsafe_poly_10 = unsafe_estimate_plot_10.Shape;
                unsafe_poly_11 = unsafe_estimate_plot_11.Shape;
                unsafe_poly_12 = unsafe_estimate_plot_12.Shape;
                unsafe_poly_13 = unsafe_estimate_plot_13.Shape;
                unsafe_poly_14 = unsafe_estimate_plot_14.Shape;
                unsafe_poly_15 = unsafe_estimate_plot_15.Shape;
                unsafe_poly_16 = unsafe_estimate_plot_16.Shape;
                unsafe_poly_17 = unsafe_estimate_plot_17.Shape;
                unsafe_poly_18 = unsafe_estimate_plot_18.Shape;
                unsafe_poly_19 = unsafe_estimate_plot_19.Shape;
                unsafe_poly_20 = unsafe_estimate_plot_20.Shape;
            end
            
            else
                if isempty(cartplot) || ~isvalid(cartplot) ...
                    || isempty(poleplot) || ~isvalid(poleplot) ...
                    || isempty(theta_lim_lt_plot)       || ~isvalid(theta_lim_lt_plot)       ...
                    || isempty(theta_lim_rt_plot)       || ~isvalid(theta_lim_rt_plot)       ...
                    || isempty(x_lim_lt_plot)           || ~isvalid(x_lim_lt_plot)           ...
                    || isempty(x_lim_rt_plot)           || ~isvalid(x_lim_rt_plot)           ...
                    || isempty(leftblockplot)           || ~isvalid(leftblockplot)           ...
                    || isempty(rightblockplot)          || ~isvalid(rightblockplot)          ...


                delete(cartplot);
                delete(poleplot);
                delete(theta_lim_lt_plot);
                delete(theta_lim_rt_plot);
                delete(x_lim_lt_plot);
                delete(x_lim_rt_plot);
                delete(leftblockplot);
                delete(rightblockplot);

                % Create polyshape objects
                cartpoly = polyshape([-cw/2,-cw/2,cw/2,cw/2],[0   ,ch       ,ch     ,0   ]);
                polepoly = polyshape([-pw/2,-pw/2,pw/2,pw/2],[ch/2,ph + ch/2,ph+ch/2,ch/2]);
                
                
                x_lim_left_block = [-x_x_lim-0.5, -x_x_lim-0.5,  -x_x_lim, -x_x_lim];
                y_lim_left_block = [-y_x_lim,      y_x_lim,      y_x_lim,  -y_x_lim ];
                x_lim_right_block = [x_x_lim+0.5,  x_x_lim+0.5,   x_x_lim,  x_x_lim];
                y_lim_right_block = [-y_x_lim,     y_x_lim,       y_x_lim, -y_x_lim];
           

                
                leftblockpoly = polyshape(x_lim_left_block,y_lim_left_block);
                rightblockpoly = polyshape(x_lim_right_block,y_lim_right_block);
                
                % Create the theta limit lines
                lineargs = {...
                    'LineStyle'     ,'--', ...
                    'Color'         ,'k' , ...
                    'LineWidth'     ,1     ...
                    };
%                 theta_lim_lt_plot = line(ha,[0 0],[0 0],lineargs{:});
%                 theta_lim_rt_plot = line(ha,[0 0],[0 0],lineargs{:});
                
                % Create the x limit lines
                x_lim_lt_plot = line(ha,[0 0],[0 0],lineargs{:});
                x_lim_rt_plot = line(ha,[0 0],[0 0],lineargs{:});
                

                cartplot = plot(ha,cartpoly,'FaceColor','b');
                cartplot.Tag = 'cartplot';
                poleplot = plot(ha,polepoly,'FaceColor','b');
                poleplot.Tag = 'poleplot';
                
                leftblockplot = plot(ha,leftblockpoly,'FaceColor','r','FaceAlpha',0.5);
                leftblockplot.Tag = 'leftblockplot';
                rightblockplot = plot(ha,rightblockpoly,'FaceColor','r','FaceAlpha',0.5);
                rightblockplot.Tag = 'rightblockplot';
               
                
                x_lim_lt_plot.Tag = 'xlimltplot';
                x_lim_rt_plot.Tag = 'xlimrtplot';
            else
                cartpoly = cartplot.Shape;
                polepoly = poleplot.Shape;
                leftblockpoly = leftblockplot.Shape;
                rightblockpoly = rightblockplot.Shape;
                
%                 safe_poly = safe_estimate_plot.Shape;
%                 unsafe_poly = unsafe_estimate_plot.Shape;
                
            end
            
            end
            
            [cartposx,~       ] = centroid(cartpoly);
            [poleposx,poleposy] = centroid(polepoly);
            
            
            dx = x - cartposx;
            dtheta = theta - atan2(cartposx-poleposx,poleposy-ch/2);
            % wrap deflection
            dtheta = atan2(sin(dtheta),cos(dtheta));
            
            cartpoly = translate(cartpoly,[dx,0]);
            polepoly = translate(polepoly,[dx,0]);
            leftblockpoly = translate(leftblockpoly,[0,0]);
            rightblockpoly = translate(rightblockpoly,[0,0]);
            
            polepoly = rotate(polepoly,rad2deg(dtheta),[x,ch/2]);
            
            cartplot.Shape = cartpoly;
            poleplot.Shape = polepoly;
            leftblockplot.Shape = leftblockpoly;
            rightblockplot.Shape = rightblockpoly;


            if traj_enable
            if n~= 0
%                 for i = 1: n/2
%                     if i<= n/6
%                         [saveposex(i),~       ] = centroid(safe_poly_set(i));
%                         safe_poly_set(i) = translate(safe_poly_set(i),[dx,0]);
% 
%                     else
%                         [unsaveposex(i-n/6),~       ] = centroid(unsafe_poly_set(i-n/6));
%                         unsafe_poly_set(i-n/6) = translate(unsafe_poly_set(i-n/6),[dx,0]);
%                     end
%                 end
                safe_poly = translate(safe_poly,[dx,0]);
                safe_poly_2 = translate(safe_poly_2,[dx,0]);
                safe_poly_3 = translate(safe_poly_3,[dx,0]);
                safe_poly_4 = translate(safe_poly_4,[dx,0]);
                safe_poly_5 = translate(safe_poly_5,[dx,0]);
                safe_poly_6 = translate(safe_poly_6,[dx,0]);
                safe_poly_7 = translate(safe_poly_7,[dx,0]);
                safe_poly_8 = translate(safe_poly_8,[dx,0]);
                safe_poly_9 = translate(safe_poly_9,[dx,0]);
                safe_poly_10 = translate(safe_poly_10,[dx,0]);
                safe_estimate_plot.Shape = safe_poly;
                safe_estimate_plot_2.Shape = safe_poly_2;
                safe_estimate_plot_3.Shape = safe_poly_3;
                safe_estimate_plot_4.Shape = safe_poly_4;
                safe_estimate_plot_5.Shape = safe_poly_5;
                safe_estimate_plot_6.Shape = safe_poly_6;
                safe_estimate_plot_7.Shape = safe_poly_7;
                safe_estimate_plot_8.Shape = safe_poly_8;
                safe_estimate_plot_9.Shape = safe_poly_9;
                safe_estimate_plot_10.Shape = safe_poly_10;

                unsafe_poly = translate(unsafe_poly,[dx,0]);
                unsafe_poly_2 = translate(unsafe_poly_2,[dx,0]);
                unsafe_poly_3 = translate(unsafe_poly_3,[dx,0]);
                unsafe_poly_4 = translate(unsafe_poly_4,[dx,0]);
                unsafe_poly_5 = translate(unsafe_poly_5,[dx,0]);
                unsafe_poly_6 = translate(unsafe_poly_6,[dx,0]);
                unsafe_poly_7 = translate(unsafe_poly_7,[dx,0]);
                unsafe_poly_8 = translate(unsafe_poly_8,[dx,0]);
                unsafe_poly_9 = translate(unsafe_poly_9,[dx,0]);
                unsafe_poly_10 = translate(unsafe_poly_10,[dx,0]);
                unsafe_poly_11 = translate(unsafe_poly_11,[dx,0]);
                unsafe_poly_12 = translate(unsafe_poly_12,[dx,0]);
                unsafe_poly_13 = translate(unsafe_poly_13,[dx,0]);
                unsafe_poly_14 = translate(unsafe_poly_14,[dx,0]);
                unsafe_poly_15 = translate(unsafe_poly_15,[dx,0]);
                unsafe_poly_16 = translate(unsafe_poly_16,[dx,0]);
                unsafe_poly_17 = translate(unsafe_poly_17,[dx,0]);
                unsafe_poly_18 = translate(unsafe_poly_18,[dx,0]);
                unsafe_poly_19 = translate(unsafe_poly_19,[dx,0]);
                unsafe_poly_20 = translate(unsafe_poly_20,[dx,0]);
                
                unsafe_estimate_plot.Shape = unsafe_poly;
                unsafe_estimate_plot_2.Shape = unsafe_poly_2;
                unsafe_estimate_plot_3.Shape = unsafe_poly_3;
                unsafe_estimate_plot_4.Shape = unsafe_poly_4;
                unsafe_estimate_plot_5.Shape = unsafe_poly_5;
                unsafe_estimate_plot_6.Shape = unsafe_poly_6;
                unsafe_estimate_plot_7.Shape = unsafe_poly_7;
                unsafe_estimate_plot_8.Shape = unsafe_poly_8;
                unsafe_estimate_plot_9.Shape = unsafe_poly_9;
                unsafe_estimate_plot_10.Shape = unsafe_poly_10;
                unsafe_estimate_plot_11.Shape = unsafe_poly_11;
                unsafe_estimate_plot_12.Shape = unsafe_poly_12;
                unsafe_estimate_plot_13.Shape = unsafe_poly_13 ;
                unsafe_estimate_plot_14.Shape = unsafe_poly_14;
                unsafe_estimate_plot_15.Shape = unsafe_poly_15;
                unsafe_estimate_plot_16.Shape = unsafe_poly_16;
                unsafe_estimate_plot_17.Shape = unsafe_poly_17;
                unsafe_estimate_plot_18.Shape = unsafe_poly_18;
                unsafe_estimate_plot_19.Shape = unsafe_poly_19;
                unsafe_estimate_plot_20.Shape = unsafe_poly_20;

            else
            end
            end
            
%             theta_lim_lt_plot.XData = [0 -x_theta_lim_] + x - pw/2;
%             theta_lim_lt_plot.YData = [0  y_theta_lim_] + ch/2;
            
%             theta_lim_rt_plot.XData = [0  x_theta_lim_] + x + pw/2;
%             theta_lim_rt_plot.YData = [0  y_theta_lim_] + ch/2;
            
            x_lim_lt_plot.XData = [0 0] - x_x_lim;
            x_lim_lt_plot.YData = [-y_x_lim y_x_lim];
            
            x_lim_rt_plot.XData = [0 0] + x_x_lim;
            x_lim_rt_plot.YData = [-y_x_lim y_x_lim];
            
            
            % if constraints are violated, change the line to red
            c1 = 'r'; c2 = [0.0 0.5 0.0];
%             if theta > theta_lim
%                 theta_lim_lt_plot.Color = c1;
%                 theta_lim_rt_plot.Color = c2;
%                 theta_lim_lt_plot.LineStyle = '-';
%                 theta_lim_rt_plot.LineStyle = '--';
%             elseif theta < -theta_lim
%                 theta_lim_lt_plot.Color = c2;
%                 theta_lim_rt_plot.Color = c1;
%                 theta_lim_lt_plot.LineStyle = '--';
%                 theta_lim_rt_plot.LineStyle = '-';
%             else
%                 theta_lim_lt_plot.Color = c2;
%                 theta_lim_rt_plot.Color = c2;
%                 theta_lim_lt_plot.LineStyle = '--';
%                 theta_lim_rt_plot.LineStyle = '--';
%             end
            if x < -x_x_lim
                x_lim_lt_plot.Color = c1;
                x_lim_rt_plot.Color = c2;
                x_lim_lt_plot.LineStyle = '-';
                x_lim_rt_plot.LineStyle = '--';
            elseif x > x_x_lim
                x_lim_lt_plot.Color = c2;
                x_lim_rt_plot.Color = c1;
                x_lim_lt_plot.LineStyle = '--';
                x_lim_rt_plot.LineStyle = '-';
            else
                x_lim_lt_plot.Color = c2;
                x_lim_rt_plot.Color = c2;
                x_lim_lt_plot.LineStyle = '--';
                x_lim_rt_plot.LineStyle = '--';
            end
            
            % Refresh rendering in figure window
            drawnow();
        end
    end
end