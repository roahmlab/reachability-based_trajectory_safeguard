%% description
close all; clear
% Note: High fidelity model state in plotted blue in dependency,
% inconsistent with the color lengend here!

% Yifei Simon Shao
% Last edited: Jan 2021
% Create and plot an example desired and braking trajectory for car each
% flag is for a different plot, toggle one of them to 1 and leave the
% others at 0
plot_time_series = 0; %plot a example of reference & failsafe trajectory, and the vehicle's simulation during tracking
plot_error_bounds =0; %plot tracking error for two different reference and two simulations. Note pink dot is centered on the reference trajectory at a time sample but the tracking error are not centered at the same around 0,0
plot_PRS = 1;% plot reference genearted by evaluating the polynomial vs reference by computing PRS zonotopes. Note the black zonotope('sliced' online), is a subset of the blue reference, which is computed in one single operation.
% Note the polynomial evaluation is a discrete representation while the
% zonotopes are overapproxiamtion, therefore ensuring safety.
%%
plot_road= true;

%% set inital condition and desired condition
x0 = 0;
y0 = 0;
phi0 =0;
vx0 =3;
delta0 = 0;
% trajectory parameters
y_des = 1;
vx_des =2;
t_peak = 2; % s
t_f = 6; % for non-braking trajectory

A_go = highway_cruising_agent() ;
% set the agents to the initial condition
z_0 = [x0;y0;phi0;vx0;delta0] ;
A_go.reset(z_0) ;

%% generate reference traj and failsafe attached at the end
[T_go,U_go,Z_go]=parameterized_traj_1d_with_brake(vx_des, vx0,phi0,y_des);
Z_go(1,:) = Z_go(1,:) + x0;%coord tf
Z_go(2,:) = Z_go(2,:) + y0;
%% Ask robot to track
A_go.move(t_f,T_go,U_go) ;

%% plotting
hi_fid_color = [0 1 0];
ref_color = [0 0 0];
set_of_ref_color = [0 1 0];
ref_set_color = [100 150 255]/255;
not_selected_ref_set_color = [0.6 0.6 0.6];
error_color= [255, 255, 180]/255;
yellow_lane_color = [255 , 187, 25]/255;

if plot_time_series
    
    h = figure('Name','Track Ref','Position',[1 1 1000 500]) ; clf;
    subplot(1,2,1)
    
    %% plot the non-braking trajectory
    axis equal ; hold on ;
    %background
    if plot_road
        fill([0 13 13 0 0],[-2.1 -2.1 6 6 -2.1],[207,207,207]/255)
        plot([0,15],[2, 2],'--w','LineWidth',8)
    end
    href = plot(Z_go(1,:),Z_go(2,:),'LineWidth',5,'Color',ref_color);
    A_go.plot();
    hhifid = plot(0,0,'LineWidth',2,'Color',hi_fid_color);
    legend([href hhifid],'Reference Trajectory','High Fidelity Trajectory');
    xlabel('x [m]')
    ylabel('y [m]')
    xlim([0 12]);
    ylim([-2 6]);
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
    
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,'CarRefTrack.pdf','-dpdf','-r0')
    
end
%% gen another trajecotry
if plot_error_bounds
    T_error = 3; % plot the error at which time instant
    
    x0 = 0;
    y0 = 0;
    phi0 =-.1;
    vx0 =2.5;
    delta0 = -0.1;
    y_des = -2;
    vx_des =3;
    t_peak = 2 ; % s
    t_f = 6;
    
    A_second = highway_cruising_agent() ;
    z_0 = [x0;y0;phi0;vx0;delta0] ;
    A_second.reset(z_0) ;
    [T_second,U_second,Z_second]=parameterized_traj_1d_with_brake(vx_des, vx0,phi0,y_des);
    
    % track the desired trajectory
    A_second.move(t_f,T_second,U_second) ;
    h2 = figure('Name','Track Ref','Position',[1 1 1300 400]) ; clf;
    subplot(1,2,1)
    %% plot
    axis equal ; hold on ;
    fill([0 14 14 0 0],[-4.1 -4.1 4 4 -4.1],[207,207,207]/255) % road color
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
    
    xlabel('x [m]')
    ylabel('y [m]')
    title('Global Coordinate')
    set(gca,'FontSize',15)
    ylim([-4,4]);xlim([0, 14])
    %%
    subplot(1,2,2);cla;hold on; grid on; axis equal
    scatter(0,0,200,[230, 195, 195]/230,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',2);
    
    [~,~,~,E_all1] = compute_2d_position_error(A_go,T_go,Z_go);
    pt_error1 = E_all1(idx_go*2+1:idx_go*2+2,:);
    % [~,idx_hi_fid]=min(abs(A_second.time - T_error));
    [~,~,~,E_all2] = compute_2d_position_error(A_second,T_second,Z_second);
    pt_error2 = E_all2(idx_second*2+1:idx_second*2+2,:);% [EloBL,EhiBL] = put_error_in_time_bins(reshape(E_all(:,1),2,[]),T_second,t_lo(:),t_hi(:));
    
    [c,center_bb,h,gen] = minBoundingBox([pt_error1 pt_error2]);
    plot(c(1,:),c(2,:),'Color',[215, 25, 28]/255,'LineWidth',3)
    error1_outer_color = [186 163 0]/255;
    error2_outer_color = [0 153 138]/255;
    scatter(pt_error1(1,:),pt_error1(2,:),80,error_color,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',2)
    scatter(pt_error2(1,:),pt_error2(2,:),80,'c','filled','MarkerEdgeColor',[0 0 0],'LineWidth',2)
    
    title('Error Coordinate')
    xlabel('x [m]')
    ylabel('y [m]')
    set(gca,'FontSize',15)
    subplot(1,2,1);
    
    scatter(pt_error1(1,:)+ Z_go(1,idx_go)       ,pt_error1(2,:)+Z_go(2,idx_go),        80,error_color,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',2)
    scatter(pt_error2(1,:)+Z_second(1,idx_second),pt_error2(2,:)+Z_second(2,idx_second),80,'c','filled','MarkerEdgeColor',[0 0 0],'LineWidth',2)
    set(h2,'Units','Inches');
    pos = get(h2,'Position');
    set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h2,'ErrorBound.pdf','-dpdf','-r0')
end

if plot_PRS
    
    y_slice_value = -0.6;% from a range of [-1,1]m, slice at a value for less convervatism
    h3 = figure('Name','Track Ref','Position',[1 1 700 500]) ; clf;
    subplot(2,1,1)
    
    %% plot the non-braking trajectory
    axis equal ; hold on ; title('Reference Trajectories')
    
    if plot_road
        fill([-3 13 13 -3 -3],[-3 -3 6 6 -3],[207,207,207]/255)
        yline(2,'-','LineWidth',7,'Color',yellow_lane_color);
        plot([-3,15],[-2, -2],'--w','LineWidth',7)
    end
    
    
    
    y_des_range = -1:0.4:1;
    
    %plot polynomial evaluation
    for y_des = y_des_range
        x0 = 0;
        y0 = 0;
        phi0 =-.1;
        vx0 = 3;
        delta0 = -0.1;
        vx_des =3;
        
        t_peak = 2 ; % s
        t_f = 6;
        A_go = highway_cruising_agent() ;
        
        z_0 = [x0;y0;phi0;vx0;delta0] ;
        A_go.reset(z_0) ;
        
        [T_go,U_go,Z_go]=parameterized_traj_1d_with_brake(vx_des, vx0,phi0,y_des);
        
        c = ref_set_color;

        if y_des ==y_slice_value
            c= ref_color;
        end
        plot(Z_go(1,:),Z_go(2,:),':','LineWidth',3,'Color',c)
    end
    
    A_go.plot(); %initail condition
    
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
    end
    set(gca,'FontSize',15);
    xlabel('x [m]');
    ylabel('y [m]');
    title('Zonotope Representation');
    
    A_second = highway_cruising_agent() ;
    
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
        
        kvc = vx0;
        kyc = y_range(y_idx);
        
        khc = phi0;
        kvd = vx_des;
        delta_c = -0.1;
        
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
        %% section 1: velocity ref goes to peak, y ref goes half way to peak
        sys = nonlinearSys(num_prs, 1, @highway_toPeak, options);
        options.x0 = [0;0;kvd;kvc;khc; kyc;0];
        %v y k_pk, kv,kh,ky, t
        options.R0 = zonotope([options.x0,diag([0;0;kvdg;kvig;h_ini_range;kyg;0])]);
        Rcont = reach(sys, options);
        
        %% section 2 velocity ref slow down, y ref goes all the way to peak
        options.R0 = Rcont{end}{1} - [0;0;0;0;0;0;2]; %initial state for reachability analysis
        options.x0 = center(options.R0);
        sys = nonlinearSys(num_prs, 1, @highway_mid, options);
        Rcont2 = reach(sys, options);
        %% section 3 velocity ref slow to 0, y ref maintain the final value
        options.R0 = Rcont2{end}{1} - [0;0;0;0;0;0;2]; %initial state for reachability analysis
        options.x0 = center(options.R0);
        sys = nonlinearSys(num_prs, 1, @highway_toStop, options);
        Rcont3 = reach(sys, options);
        
        Rcont_all = [Rcont;Rcont2;Rcont3];
        
        for i = 1:2:length(Rcont_all)
            %skip a few for clarity
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



