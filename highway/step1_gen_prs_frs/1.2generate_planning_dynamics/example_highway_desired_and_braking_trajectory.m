%% description
close all; clear 
% Create and plot an example desired and braking trajectory for the segway
plot_agent = 0;
plot_time_series = 0;
plot_error_bounds =0 ;
plot_PRS = 1;
plot_online= 0;
% Author: Shreyas Kousik
% Created: 9 Mar 2020
% Updated: 10 Mar 2020
%
%% user parameters
% robot initial condition
% w_0 = 0.0 ; % rad/s
% v_0 = 0.75 ; % m/s
plot_road= true;

x0 = 0;
y0 = 0;
phi0 =0;
vx0 =3;
delta0 = 0;
% trajectory parameters
y_des = 1;
vx_des =2;
% w_des = 1.0 ; % rad/s
% v_des = 0.25 ; % m/s
t_peak = 2; % s
t_f = 6; % for non-braking trajectory

%% automated from here
% make segway for each trajectory
A_go = highway_cruising_agent() ;
% need to add this to agent constructor to work properly
A.footprint_vertices = [-2.4,-1.5,-1.5 0 0.3     2    2   2.4 2.4 2  2   0.3 0 -1.5 -1.5 -2.4 -2.4;
                -0.5,-0.5 -1   -1 -0.5    -0.3 -1   -1  1  1 0.3  0.5 1 1   0.5 0.5 -0.5];
% A.plot_footprint_color = [1 0 0];
% A_brk = segway_agent() ;

% set the agents to the initial condition
z_0 = [x0;y0;phi0;vx0;delta0] ;
A_go.reset(z_0) ;
% A_brk.reset(z_0) ;

% get stopping time
% t_stop = A_go.max_speed ./ A_go.max_accel ;

% make desired trajectory
% [T_go,U_go,Z_go] = make_highway_desired_trajectory(t_f,y_des,vx_des,z_0) ;
% [T_go,U_go,Z_go] = make_segway_desired_trajectory(t_f,w_des,v_des) ;

[T_go,U_go,Z_go]=parameterized_traj_1d_with_brake(vx_des, vx0,phi0,y_des);
Z_go(1,:) = Z_go(1,:) + x0;
Z_go(2,:) = Z_go(2,:) + y0;

% make braking trajectory
% [T_brk,U_brk,Z_brk] = make_segway_braking_trajectory(t_plan,t_stop,w_des,v_des) ;

% track the desired trajectory
A_go.move(t_f,T_go,U_go) ;
% Z_go = zeros()
% A_brk.move(t_plan+t_stop,T_brk,U_brk,Z_brk) ;

%% plotting
hi_fid_color = [0 1 0];
ref_color = [0 0 0];
set_of_ref_color = [0 1 0];
ref_set_color = [100 150 255]/255;
not_selected_ref_set_color = [0.6 0.6 0.6];
error_color= [255, 255, 180]/255;
yellow_lane_color = [255 , 187, 25]/255;
% error_color= [102, 178, 180]/255;


if plot_agent
    
    x0 = 0;
y0 = 0;
phi0 =0;
vx0 =0;
delta0 = 0;
% trajectory parameters
y_des = 0;
vx_des =0;
% w_des = 1.0 ; % rad/s
% v_des = 0.25 ; % m/s
t_peak = 2; % s
t_f = 6; % for non-braking trajectory

%% automated from here
% make segway for each trajectory
A_go = highway_cruising_agent() ;
% A.plot_footprint_color = [1 0 0];
% A_brk = segway_agent() ;

% set the agents to the initial condition
z_0 = [x0;y0;phi0;vx0;delta0] ;
A_go.reset(z_0) ;
% A_brk.reset(z_0) ;

% get stopping time
% t_stop = A_go.max_speed ./ A_go.max_accel ;

% make desired trajectory
% [T_go,U_go,Z_go] = make_highway_desired_trajectory(t_f,y_des,vx_des,z_0) ;
% [T_go,U_go,Z_go] = make_segway_desired_trajectory(t_f,w_des,v_des) ;

[T_go,U_go,Z_go]=parameterized_traj_1d_with_brake(vx_des, vx0,phi0,y_des);


% make braking trajectory
% [T_brk,U_brk,Z_brk] = make_segway_braking_trajectory(t_plan,t_stop,w_des,v_des) ;

% track the desired trajectory
A_go.move(t_f,T_go,U_go) ;

h = figure('Name','Agent','Position',[1 1 2000 600]) ;
 axis equal ; hold on ; 
 A_go.plot();
 
end

if plot_time_series
    
    h = figure('Name','Track Ref','Position',[1 1 1000 500]) ; clf;
    % iax = 1; % Or whichever
    % subaxis(4, 6, iax, 'sh', 0.03, 'sv', 0.01, 'padding', 0, 'margin', 0);
     subplot(1,2,1)
    
    %% plot the non-braking trajectory
    axis equal ; hold on ; 
    %background
    if plot_road 
        fill([0 13 13 0 0],[-2.1 -2.1 6 6 -2.1],[207,207,207]/255)
%         yline(-2,'-','LineWidth',8,'Color',[247, 207, 25]/255);
        plot([0,15],[2, 2],'--w','LineWidth',8)
%         yline(2,'--w','LineWidth',10);
    end
%     set(gca,'Color',);
    href = plot(Z_go(1,:),Z_go(2,:),'LineWidth',5,'Color',ref_color);
    A_go.plot();
    hhifid = plot(0,0,'LineWidth',2,'Color',hi_fid_color);
    %plot obstacle
%  
    % legend(hA(3),'hi fid');
    legend([href hhifid],'Reference Trajectory','High Fidelity Trajectory');
    % plot_path(Z_go(1:2,:),'b--','linewidth',1.5)
    xlabel('x [m]')
    ylabel('y [m]')
    xlim([0 12]);
    ylim([-2 6]);
    % title('Workspace Reference and Tracking Trajectory')
    set(gca,'FontSize',15)
    
    subplot(2,2,2);cla;hold on; grid on;
    plot(T_go,U_go(1,:),'LineWidth',5,'Color',ref_color)
    plot(A_go.time,A_go.state(4,:),'LineWidth',2, 'Color',hi_fid_color);
    ylabel('Veclocity [m/s]')
    set(gca,'FontSize',15)
    
    subplot(2,2,4);cla;hold on; grid on;
    plot(T_go,U_go(2,:),'LineWidth',5,'Color',ref_color)%
    plot(A_go.time,A_go.state(2,:),'LineWidth',2,'Color',hi_fid_color);
    xlabel('Time [s]')
    ylabel('y Position [m]')
    set(gca,'FontSize',15)
%     set(gcf, 'InvertHardcopy', 'off')
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,'CarRefTrack.pdf','-dpdf','-r0')
    
end
%% gen another trajecotry
if plot_error_bounds
    T_error = 3;
    
    x0 = 0;
    y0 = 0;
    phi0 =-.1;
    vx0 =2.5;
    delta0 = -0.1;
    % trajectory parameters
    y_des = -2;
    vx_des =3;
    % w_des = 1.0 ; % rad/s
    % v_des = 0.25 ; % m/s
    t_peak = 2 ; % s
    t_f = 6; % for non-braking trajectory
    
    %% automated from here
    % make segway for each trajectory
    A_second = highway_cruising_agent() ;
    % A.plot_footprint_color = [1 0 0];
    % A_brk = segway_agent() ;
    
    % set the agents to the initial condition
    z_0 = [x0;y0;phi0;vx0;delta0] ;
    A_second.reset(z_0) ;
    % A_brk.reset(z_0) ;
    
    % get stopping time
    % t_stop = A_go.max_speed ./ A_go.max_accel ;
    
    % make desired trajectory
    % [T_go,U_go,Z_go] = make_highway_desired_trajectory(t_f,y_des,vx_des,z_0) ;
    % [T_go,U_go,Z_go] = make_segway_desired_trajectory(t_f,w_des,v_des) ;
    
    [T_second,U_second,Z_second]=parameterized_traj_1d_with_brake(vx_des, vx0,phi0,y_des);
    
    
    
    % make braking trajectory
    % [T_brk,U_brk,Z_brk] = make_segway_braking_trajectory(t_plan,t_stop,w_des,v_des) ;
    
    % track the desired trajectory
    A_second.move(t_f,T_second,U_second) ;
    h2 = figure('Name','Track Ref','Position',[1 1 1300 400]) ; clf;
    % iax = 1; % Or whichever
    % subaxis(4, 6, iax, 'sh', 0.03, 'sv', 0.01, 'padding', 0, 'margin', 0);
    subplot(1,2,1)
    %% plot the non-braking trajectory
    axis equal ; hold on ;
    fill([0 14 14 0 0],[-4.1 -4.1 4 4 -4.1],[207,207,207]/255) % road color
%      set(gca,'Color',[207,207,207]/255);
%     line1 = yline(4,'-','LineWidth',7,'Color',yellow_lane_color); % lane
    yline(0,'--w','LineWidth',7);
    href = plot(Z_go(1,:),Z_go(2,:),'LineWidth',5,'Color',ref_color); %reference
    [~,idx_go] = min(abs(T_go - T_error));
    A_go.plot_at_time(T_error); % agent and high fid
    hhifid = plot(0,0,'LineWidth',2,'Color',hi_fid_color);
    ref_point_outer_color = [221 127 127]/255;
    scatter(Z_go(1,idx_go),Z_go(2,idx_go),200,[230, 195, 195]/230,'filled','MarkerEdgeColor',[0  0  0],'LineWidth',2)
    %ref point
    
    
    plot(Z_second(1,:),Z_second(2,:),'LineWidth',5,'Color',ref_color) 
    A_second.plot_at_time(T_error);
    hhifid = plot(0,0,'LineWidth',2,'Color',hi_fid_color);
    [~,idx_second] = min(abs(T_second - T_error));
    scatter(Z_second(1,idx_second),Z_second(2,idx_second),200,[230, 195, 195]/230,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',2)
    
    % plot_path(Z_go(1:2,:),'b--','linewidth',1.5)
    xlabel('x [m]')
    ylabel('y [m]')
    title('Global Coordinate')
    % title('Workspace Reference and Tracking Trajectory')
    set(gca,'FontSize',15)
    ylim([-4,4]);xlim([0, 14])
    %%
    subplot(1,2,2);cla;hold on; grid on; axis equal
    scatter(0,0,200,[230, 195, 195]/230,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',2);
    % t_lo= 0:0.1;
    % t_hi= dt:dt:t_total;
    % [~,idx_hi_fid]=min(abs(A_go.time - T_error));
    [~,~,~,E_all1] = compute_2d_position_error(A_go,T_go,Z_go);
    pt_error1 = E_all1(idx_go*2+1:idx_go*2+2,:);
    % [~,idx_hi_fid]=min(abs(A_second.time - T_error));
    [~,~,~,E_all2] = compute_2d_position_error(A_second,T_second,Z_second);
    pt_error2 = E_all2(idx_second*2+1:idx_second*2+2,:);% [EloBL,EhiBL] = put_error_in_time_bins(reshape(E_all(:,1),2,[]),T_second,t_lo(:),t_hi(:));
    % [EloBR,EhiBR] = put_error_in_time_bins(reshape(E_all(:,2),2,[]),T_second,t_lo(:),t_hi(:));
    % [EloTR,EhiTR] = put_error_in_time_bins(reshape(E_all(:,3),2,[]),T_second,t_lo(:),t_hi(:));
    % [EloTL,EhiTL] = put_error_in_time_bins(reshape(E_all(:,4),2,[]),T_second,t_lo(:),t_hi(:));
    % e_verts = [EloBL EhiBL EloBR EhiBR EloTR EhiTR EloTL EhiTL];
    % verts=reshape(error_verts(i,:)',2,[]);
    % [c,center_bb,h,gen] = minBoundingBox(verts);
    
    [c,center_bb,h,gen] = minBoundingBox([pt_error1 pt_error2]);
    plot(c(1,:),c(2,:),'Color',[215, 25, 28]/255,'LineWidth',3)
    error1_outer_color = [186 163 0]/255;
    error2_outer_color = [0 153 138]/255;
    scatter(pt_error1(1,:),pt_error1(2,:),80,error_color,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',2)
    scatter(pt_error2(1,:),pt_error2(2,:),80,'c','filled','MarkerEdgeColor',[0 0 0],'LineWidth',2)
    
    title('Error Coordinate')
    xlabel('x [m]')
    ylabel('y [m]')
%     scatter(Z_go(1,idx_go),Z_go(2,idx_go),200,[230, 195, 195]/230,'filled','MarkerEdgeColor',[0  0  0],'LineWidth',2);
    set(gca,'FontSize',15)
    subplot(1,2,1);
    
    scatter(pt_error1(1,:)+ Z_go(1,idx_go)       ,pt_error1(2,:)+Z_go(2,idx_go),        80,error_color,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',2)
    scatter(pt_error2(1,:)+Z_second(1,idx_second),pt_error2(2,:)+Z_second(2,idx_second),80,'c','filled','MarkerEdgeColor',[0 0 0],'LineWidth',2)
%     set(gcf, 'InvertHardcopy', 'off')
    set(h2,'Units','Inches');
    pos = get(h2,'Position');
    set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h2,'ErrorBound.pdf','-dpdf','-r0')
end

if plot_PRS
    
    y_slice_value = -0.6;
    h3 = figure('Name','Track Ref','Position',[1 1 700 500]) ; clf;
    subplot(2,1,1)
    
    %% plot the non-braking trajectory
    axis equal ; hold on ;   title('Reference Trajectories')

    if plot_road 
        fill([-3 13 13 -3 -3],[-3 -3 6 6 -3],[207,207,207]/255)
        yline(2,'-','LineWidth',7,'Color',yellow_lane_color);
        plot([-3,15],[-2, -2],'--w','LineWidth',7)
%         yline(2,'--w','LineWidth',10);
    end
    for i = 2
        if i == 1
            y_des_range = -3:0.4:-1;
        elseif i == 2
            y_des_range = -1:0.4:1;
        elseif i == 3
            y_des_range = 1:0.4:3;
        end
        
        for y_des = y_des_range
            x0 = 0;
            y0 = 0;
            phi0 =-.1;
            vx0 = 3;
            delta0 = -0.1;
            % trajectory parameters
            %         y_des = -2;
            vx_des =3;
            % w_des = 1.0 ; % rad/s
            % v_des = 0.25 ; % m/s
            t_peak = 2 ; % s
            t_f = 6; % for non-braking trajectory
            A_go = highway_cruising_agent() ;
            % A.plot_footprint_color = [1 0 0];
            % A_brk = segway_agent() ;
            
            % set the agents to the initial condition
            z_0 = [x0;y0;phi0;vx0;delta0] ;
            A_go.reset(z_0) ;
            % A_brk.reset(z_0) ;
            
            % get stopping time
            % t_stop = A_go.max_speed ./ A_go.max_accel ;
            
            % make desired trajectory
            % [T_go,U_go,Z_go] = make_highway_desired_trajectory(t_f,y_des,vx_des,z_0) ;
            % [T_go,U_go,Z_go] = make_segway_desired_trajectory(t_f,w_des,v_des) ;
            
            [T_go,U_go,Z_go]=parameterized_traj_1d_with_brake(vx_des, vx0,phi0,y_des);
            if i == 2
                c = ref_set_color;
            else
                c = not_selected_ref_set_color;
            end
            if y_des ==y_slice_value
                c= ref_color;
            end
            plot(Z_go(1,:),Z_go(2,:),':','LineWidth',3,'Color',c)
        end
    end
    % iax = 1; % Or whichever
    % subaxis(4, 6, iax, 'sh', 0.03, 'sv', 0.01, 'padding', 0, 'margin', 0);
    
    %     href = plot(Z_go(1,:),Z_go(2,:),'LineWidth',5,'Color',ref_color)
    A_go.plot();
    %     hhifid = plot(0,0,'LineWidth',2,'Color',hi_fid_color)
    % legend(hA(3),'hi fid');
    %     legend([href hhifid],'Reference Trajectory','High Fidelity Trajectory');
    % plot_path(Z_go(1:2,:),'b--','linewidth',1.5)
    xlabel('x [m]')
    ylabel('y [m]')
    % title('Workspace Reference and Tracking Trajectory')
    set(gca,'FontSize',15)
    xlim([-3 13]);ylim([-2.3 2])
   
    
    %%
    subplot(2,1,2); hold on; axis equal; 
    if plot_road 
        fill([-3 13 13 -3 -3],[-3 -3 6 6 -3],[207,207,207]/255)
        yline(2,'-','LineWidth',7,'Color',yellow_lane_color);
        plot([-3,15],[-2, -2],'--w','LineWidth',7)
%         yline(2,'--w','LineWidth',10);
    end
    set(gca,'FontSize',15);
    xlabel('x [m]');
    ylabel('y [m]');
    title('Zonotope Representation');
    
    A_second = highway_cruising_agent() ;
    % A.plot_footprint_color = [1 0 0];
    % A_brk = segway_agent() ;
    
    % set the agents to the initial condition
    z_0 = [x0;y0;phi0;vx0;delta0] ;
    A_second.reset(z_0) ;
    
    A_second.plot();
    
    dt = 0.05;
    dt_stop = 0.1;
    t_peak = 2;
    t_mid = 2;
    t_end = 2;
    t_total = 6;
    
    y_range = [-2 0 2];
    del_ini_range = 0.1; %0.07; %from  experiment
    h_ini_range = 0.05;
    kyg = 1;
    kvdg = 1;
    kvig = 1;
    
    for y_idx = [1 3 2]
%                     my_counter= my_counter +1
                    
                    
        kvc = vx0;
        kyc = y_range(y_idx);

        khc = phi0;
        kvd = vx_des;
        delta_c = -0.1;
        %                     kvc = 6.5;
        %                     kvd = 7;
        %                     figure(1);hold on
        %                         clf;hold on;
        % for i = 1:length(Rcont)
        %     p2 = plotFilled(zonotope_slice(Rcont{i}{1}, slice_dim, slice_parameter), [3, 4], 'b');
        %     p2.FaceAlpha = 0.25;
        % end
        % plot(y(:,3),y(:,4),'r');
        %% section 1 prs
        num_prs = 7;
        options = struct;
        options.tStart = 0;
        options.tFinal = t_peak;
        options.timeStep=dt; %time step size for reachable set computation
        options.taylorTerms=3; %number of taylor terms for reachable sets
        options.zonotopeOrder= 3; %zonotope order... increase this for more complicated systems.
        options.maxError =9e100*ones(num_prs, 1); % this controls splitting, set it high to avoid splitting
        options.verbose = 0;

        options.uTrans = 0; % set input = 0
        options.U = zonotope([0,0 ]);
        options.advancedLinErrorComp = 0;
        options.tensorOrder = 1;
        options.reductionInterval = inf;
        options.reductionTechnique = 'girard';

        sys = nonlinearSys(num_prs, 1, @highway_toPeak, options);
        options.x0 = [0;0;kvd;kvc;khc; kyc;0];
        %v y k_pk, kv,kh,ky, t
        options.R0 = zonotope([options.x0,diag([0;0;kvdg;kvig;h_ini_range;kyg;0])]);
        Rcont = reach(sys, options);

        %% section 2 prs
        %                     options.tStart = 0; %start time
        %                     options.tFinal = 2; %final time
        % parameters: [p; v0; a0; v_peak; t];
        %                     p_peak = Rcont{end}{1}.Z(1, 1);
        options.R0 = Rcont{end}{1} - [0;0;0;0;0;0;2]; %initial state for reachability analysis
        options.x0 = center(options.R0);
        sys = nonlinearSys(num_prs, 1, @highway_mid, options);
        Rcont2 = reach(sys, options);
        %% section 3 stop
        options.R0 = Rcont2{end}{1} - [0;0;0;0;0;0;2]; %initial state for reachability analysis
        options.x0 = center(options.R0);
        sys = nonlinearSys(num_prs, 1, @highway_toStop, options);
        Rcont3 = reach(sys, options);
        
        Rcont_all = [Rcont;Rcont2;Rcont3];
        
        for i = 1:2:length(Rcont_all)
%             if i > 80 && mod(i-1,8)
%                 continue;
%             end
            if y_idx == 2
            prs = project(zonotope_slice(Rcont_all{i}{1}, [3;4;5], [vx_des;vx0;phi0]), [1, 2]);
            p3= plotFilled(prs,[1,2]);
            p3.FaceAlpha = 0.3;
            
                p3.FaceColor = ref_set_color;
                p3.EdgeColor = ref_set_color;
                prs_sliced = project(zonotope_slice(Rcont_all{i}{1}, [3;4;5;6], [vx_des;vx0;phi0;y_slice_value]), [1, 2]);
                p4= plotFilled(prs_sliced,[1,2]);
                p4.FaceAlpha = 0.3;
                p4.FaceColor = ref_color;
                p4.EdgeColor = ref_color;
            else
                p3.FaceColor = not_selected_ref_set_color;
                p3.EdgeColor = not_selected_ref_set_color;
            end
        end
    end
    
 xlim([-3 13]);ylim([-2.3 2])
% set(gcf, 'InvertHardcopy', 'off')
set(h3,'Units','Inches');
pos = get(h3,'Position');
set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h3,'RefvsZono.pdf','-dpdf','-r0')    
    
end


% if plot_online
%     T_error = 0;
%      O = make_box([4,2]);
%     O = O + repmat([13;2],1,5);
%     O2 = make_box([15 0.5]);
%      O2 = O2 + repmat([15/2;4],1,5);
%     % A.plot_footprint_color = [1 0 0];
%     % A_brk = segway_agent() ;
%     
%     % set the agents to the initial condition
%    
%     % A_brk.reset(z_0) ;
%     
%     % get stopping time
%     % t_stop = A_go.max_speed ./ A_go.max_accel ;
%     
%     % make desired trajectory
%     % [T_go,U_go,Z_go] = make_highway_desired_trajectory(t_f,y_des,vx_des,z_0) ;
%     % [T_go,U_go,Z_go] = make_segway_desired_trajectory(t_f,w_des,v_des) ;
%     
% %     [T_second,U_second,Z_second]=parameterized_traj_1d_with_brake(vx_des, vx0,phi0,y_des);
%     
%     
%     
%     % make braking trajectory
%     % [T_brk,U_brk,Z_brk] = make_segway_braking_trajectory(t_plan,t_stop,w_des,v_des) ;
%     
%     % track the desired trajectory
% %     A_second.move(t_f,T_second,U_second) ;
%     h2 = figure('Name','Track Ref','Position',[1 1 1300 400]) ; clf;
%     
%     % iax = 1; % Or whichever
%     % subaxis(4, 6, iax, 'sh', 0.03, 'sv', 0.01, 'padding', 0, 'margin', 0);
%     subplot(1,2,1)
%     
%     %% plot the non-braking trajectory
%     axis equal ; hold on ;
%     
%     fill([0 14 14 0 0],[-4.1 -4.1 4 4 -4.1],[207,207,207]/255)
%     hobs = fill(O(1,:),O(2,:),[1 0 0]);
%      set(hobs,'facealpha',.5)
%      hobs = fill(O2(1,:),O2(2,:),[1 0 0]);
%      set(hobs,'facealpha',.5)
% %      set(gca,'Color',[207,207,207]/255);
%     line1 = yline(-4,'-','LineWidth',10,'Color',[247, 207, 25]/255);
%     yline(0,'--w','LineWidth',10);
%     href = plot(Z_go(1,:),Z_go(2,:),'LineWidth',5,'Color',ref_color);
%     [~,idx_go] = min(abs(T_go - T_error));
%     A_go.plot_at_time(T_error);
%     hhifid = plot(0,0,'LineWidth',2,'Color',hi_fid_color);
% %     scatter(Z_go(1,idx_go),Z_go(2,idx_go),200,[230, 195, 195]/230,'filled','MarkerEdgeColor',[221 127 127]/255,'LineWidth',2)
%     
%     % legend(hA(3),'hi fid');
%     % legend([href hhifid],'Reference Trajectory','High Fidelity Trajectory');
%     
% %     plot(Z_second(1,:),Z_second(2,:),'LineWidth',5,'Color',ref_color)
% %     A_second.plot_at_time(T_error);
%     hhifid = plot(0,0,'LineWidth',2,'Color',hi_fid_color);
% %     [~,idx_second] = min(abs(T_second - T_error));
% %     scatter(Z_second(1,idx_second),Z_second(2,idx_second),200,[230, 195, 195]/230,'filled','MarkerEdgeColor',[221 127 127]/255,'LineWidth',2)
%     
%     % plot_path(Z_go(1:2,:),'b--','linewidth',1.5)
%     xlabel('x [m]')
%     ylabel('y [m]')
%     title('Global Coordinate')
%     % title('Workspace Reference and Tracking Trajectory')
%     set(gca,'FontSize',15)
%     ylim([-4,4]);xlim([0, 14])
%     %%
%     subplot(1,2,2);cla;hold on; grid on; axis equal
% %     scatter(0,0,200,[230, 195, 195]/230,'filled','MarkerEdgeColor',[221 127 127]/255,'LineWidth',2);
%     % t_lo= 0:0.1;
%     % t_hi= dt:dt:t_total;
%     % [~,idx_hi_fid]=min(abs(A_go.time - T_error));
%     [~,~,~,E_all1] = compute_2d_position_error(A_go,T_go,Z_go);
%     pt_error1 = E_all1(idx_go*2+1:idx_go*2+2,:);
%     % [~,idx_hi_fid]=min(abs(A_second.time - T_error));
% %     [~,~,~,E_all2] = compute_2d_position_error(A_second,T_second,Z_second);
% %     pt_error2 = E_all2(idx_second*2+1:idx_second*2+2,:);% [EloBL,EhiBL] = put_error_in_time_bins(reshape(E_all(:,1),2,[]),T_second,t_lo(:),t_hi(:));
%     % [EloBR,EhiBR] = put_error_in_time_bins(reshape(E_all(:,2),2,[]),T_second,t_lo(:),t_hi(:));
%     % [EloTR,EhiTR] = put_error_in_time_bins(reshape(E_all(:,3),2,[]),T_second,t_lo(:),t_hi(:));
%     % [EloTL,EhiTL] = put_error_in_time_bins(reshape(E_all(:,4),2,[]),T_second,t_lo(:),t_hi(:));
%     % e_verts = [EloBL EhiBL EloBR EhiBR EloTR EhiTR EloTL EhiTL];
%     % verts=reshape(error_verts(i,:)',2,[]);
%     % [c,center_bb,h,gen] = minBoundingBox(verts);
%     
% %     [c,center_bb,h,gen] = minBoundingBox([pt_error1 pt_error2]);
% %     plot(c(1,:),c(2,:),'Color',[215, 25, 28]/255,'LineWidth',3)
%     scatter(pt_error1(1,:),pt_error1(2,:),80,error_color,'filled','MarkerEdgeColor',[186 163 0]/255,'LineWidth',2)
%     scatter(pt_error2(1,:),pt_error2(2,:),80,'c','filled','MarkerEdgeColor',[0 153 138]/255,'LineWidth',2)
%     
%     title('Error Coordinate')
%     xlabel('x [m]')
%     ylabel('y [m]')
%     set(gca,'FontSize',15)
%     subplot(1,2,1);
%     
%     scatter(pt_error1(1,:)+ Z_go(1,idx_go)       ,pt_error1(2,:)+Z_go(2,idx_go),        80,error_color,'filled','MarkerEdgeColor',[186 163 0]/255,'LineWidth',2)
%     scatter(pt_error2(1,:)+Z_second(1,idx_second),pt_error2(2,:)+Z_second(2,idx_second),80,'c','filled','MarkerEdgeColor',[0 153 138]/255,'LineWidth',2)
% %     set(gcf, 'InvertHardcopy', 'off')
%     set(h2,'Units','Inches');
%     pos = get(h2,'Position');
%     set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     print(h2,'OnlineOps.pdf','-dpdf','-r0')
% end

