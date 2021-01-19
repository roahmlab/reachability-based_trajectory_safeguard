% In 1.1, 1.2, we defined the high fidelity model and the planning model, and a reference version and a CORA version of the planning model
% Here we use CORA to compute PRS and simulation to compute ERS. 
clear; clc;close all
plot_flag = 1; % Plot Error for samples in a time bin and bound it using a rectangle to save space, then plot a simulation and sliced FRS in red.
debug_flag = 0; % debug flag plot the PRS sliced at initial condition, and plot the reference in blue, and simulation trajecotry in red.
w = warning ('off','all');

v_range =2:2:4;% set the centers of the velocity, y_desired, heading and delta
%v = 1---3 , 3----5
v_des_range = v_range;
y_range = [0] ; % +-1 can split into more partitions to increase accuracy
h_range = -0.3:0.1:0.3; % this range is from simulation of random expriments. needs to encompass all initial condition sets
del_range = [-0.2 0  0.2 ];

%set the size of the generators
del_ini_range = 0.1; %wheel angle uncertainty
h_ini_range = 0.05; %heading uncertainty
kyg = 1; %desired y uncertainty
kvdg = 1; %desired v can be different for different speeds, keep it same for now
kvig = 1;%initial v uncertainty

num_v_des =length(v_range);

yr = 1:length(y_range); % set sizes for loop
vr = 1:length(v_range);
hr = 1:length(h_range);
deltar = 1:length(del_range);

dt = 0.05; %dt for zonotope time samples
t_peak = 2; % 6 seconds of total time for each plan, 2 reaching a desired velocity, then reducing speed to 0 after
t_mid = 2; % 4 seconds to reach the desired y position, y_position held constant after
t_end = 2; 
t_total = 6;
t_lo= 0:dt:t_total-dt; % time sections to contain error data
t_hi= dt:dt:t_total;
N_t_bins = length(t_lo);



tic
%%
res=cell(length(v_range),length(y_range),length(h_range),length(del_range),num_v_des,3); % each section has 2 seconds of zonotopes, 3 sections in total
error_table = cell(length(v_range),length(y_range),length(h_range),length(del_range),num_v_des);

A = highway_cruising_agent;
A.integrator_time_discretization = 0.03;
A.integrator_type = 'ode45';
if plot_flag
    figure('Renderer', 'painters', 'Position', [10 10 800*1.5 400*1.5])
end

vd_idx_max_a = zeros(length(v_range),1);
vd_idx_min_a = zeros(length(v_range),1);
pts = [];

% set limits on the range of avaliable commanded velocity, current velocity - 4 to
% v_cur + 2
for v_idx = vr
    vi_outer = v_range(v_idx);
    vd_idx_max = find(v_range==vi_outer+2);
    if isempty(vd_idx_max)
        vd_idx_max = length(v_range);
    end
    vd_idx_min = find(v_range==vi_outer-4);
    if isempty(vd_idx_min)
        vd_idx_min = 1;
    end
    vd_idx_max_a(v_idx)=vd_idx_max;
    vd_idx_min_a(v_idx)=vd_idx_min;
    for vd_idx_plot=vd_idx_min:vd_idx_max
        pts = [pts [v_idx; vd_idx_plot]];
    end
end


my_counter = 0;
rng('shuffle')
% for each range of initial condition
for v_idx = vr 
    toc
    v_range(v_idx)
    for h_idx =hr
        for vd_idx =vd_idx_min_a(v_idx) :vd_idx_max_a(v_idx)
            for del_idx = deltar
                for y_idx = yr
                    % for each range of initial condition, generate PRS
                    my_counter= my_counter +1
                    
                    
                    kvc = v_range(v_idx); % set zonotope centers
                    kyc = y_range(y_idx);
                    
                    khc = h_range(h_idx);
                    kvd = v_range(vd_idx);
                    delta_c = del_range(del_idx);
                   
                    %% section 1 to peak
                             % states| parameters   | time
                    num_prs = 7; %v y k_pk, kv,kh,ky, t   
                    %initial conitions are treated the same as target
                    %parameters. During online when slicing, initial
                    %conditions are directly sliced with initial conditions, but target ones are
                    %evaluated at samples.
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
                                 %v y v_desired@2sec, vini, hini, y_desired@4sec, t
                    options.R0 = zonotope([options.x0,diag([0;0;kvdg;kvig;h_ini_range;kyg;0])]);
                    Rcont = reach(sys, options);
                    
                    %% section 2 mid 
                    % set final condition of previous zonotope as initial 
                    % condition of this one, also decrease time
                    options.R0 = Rcont{end}{1} - [0;0;0;0;0;0;2]; %initial state for reachability analysis
                    options.x0 = center(options.R0);
                    sys = nonlinearSys(num_prs, 1, @highway_mid, options);
                    Rcont2 = reach(sys, options);
                    
                    %% section 3 stop
                    options.R0 = Rcont2{end}{1} - [0;0;0;0;0;0;2]; %initial state for reachability analysis
                    options.x0 = center(options.R0);
                    sys = nonlinearSys(num_prs, 1, @highway_toStop, options);
                    Rcont3 = reach(sys, options);
                    
                    %% gen ERS
                    % do simulation on edge of the each parameter and
                    % calculate error. PRS overapproximate ref traj,
                    % ERS overapproximate tracking error wrt ref traj, so
                    % their sum overapproximate FRS
                    
                    Ex_lo_cur_vert = inf*zeros(N_t_bins,1) ;
                    Ey_lo_cur_vert = inf*zeros(N_t_bins,1) ;
                    %     Ez_lo_cur = zeros(N_t_bins,1) ;
                    Ex_hi_cur_vert = -inf*zeros(N_t_bins,1) ;
                    Ey_hi_cur_vert = -inf*zeros(N_t_bins,1) ;
                    
                    Ex_lo_cur = inf*zeros(N_t_bins,1) ;
                    Ey_lo_cur = inf*zeros(N_t_bins,1) ;
                    %     Ez_lo_cur = zeros(N_t_bins,1) ;
                    Ex_hi_cur = -inf*zeros(N_t_bins,1) ;
                    Ey_hi_cur = -inf*zeros(N_t_bins,1) ;
                    
                    error_verts = zeros(N_t_bins,0);
                    error_data_cur = zeros(5,N_t_bins);
                    
                    %% sim high fidelity to get error
                    num_sims = 3; eps =   0.01;
                    % simulate 3 sample points inside the range of
                    % parameters.
                    for kvc_sim = linspace(kvc-kvig+eps, kvc+kvig-eps, num_sims)
                        for kyc_sim = linspace(kyc-kyg+eps, kyc+kyg-eps, num_sims)
                            for khc_sim = linspace(khc-h_ini_range+eps,khc+h_ini_range-eps,num_sims)
                                for kvd_sim = linspace(kvd - kvdg+eps,kvd + kvdg-eps,num_sims )
                                    for delta_c_sim = linspace(delta_c - del_ini_range+eps,delta_c + del_ini_range-eps,num_sims )
                                        ini_cond = [0;0;khc_sim;kvc_sim;delta_c_sim];
                                        A.reset(ini_cond);
                                        [T_ref,U_ref,Z_ref]=parameterized_traj_1d_with_brake(kvd_sim, kvc_sim,khc_sim,kyc_sim);
                                        %T = time; U_ref = reference for
                                        %LLC, vx and y; Z_ref = position
                                        %reference for x and y. to
                                        %calculate error in x and y dim
                                        A.move(t_total,T_ref,U_ref,Z_ref) ;
                                        if debug_flag
                                            A.plot();
                                            clf;hold on;axis auto;plot(Z_ref(1,:),Z_ref(2,:));plot(A.state(1,:),A.state(2,:))
                                            for i = 1:length(Rcont)
                                                z = project(zonotope_slice(Rcont{i}{1}, [3;4;5;6], [kvd_sim;kvc_sim;khc_sim;kyc_sim]), [1, 2]);
                                                p3= plot(z,[1,2],'r');
                                            end
                                            peak_count = length(Rcont);
                                            for i = 1:length(Rcont2)
                                                z = project(zonotope_slice(Rcont2{i}{1}, [3;4;5;6], [kvd_sim;kvc_sim;khc_sim;kyc_sim]), [1, 2]);
                                                p3= plot(z,[1,2],'g');%
                                            end
                                            two_sec_count = peak_count + length(Rcont2);
                                            for i = 1:length(Rcont3)
                                                %  [0;0;kvdg;kvig;h_ini_range;kyg;0]
                                                z = project(zonotope_slice(Rcont3{i}{1}, [3;4;5;6], [kvd_sim;kvc_sim;khc_sim;kyc_sim]), [1, 2]);
                                                p3= plot(z,[1,2],'b');
                                                %                             p3.FaceAlpha = 0.05;
                                        
                                            end
                                            drawnow
                                        end
                                        %                                         plot_arr =cell(1);
                                        [E_cur,Ehi,Elo,E_all] = compute_2d_position_error(A,T_ref,Z_ref);
                                        % get error at all reference points
                                        % and place them in time bins
                                        [EloBL,EhiBL] = put_error_in_time_bins_2d(reshape(E_all(:,1),2,[]),T_ref,t_lo(:),t_hi(:));
                                        [EloBR,EhiBR] = put_error_in_time_bins_2d(reshape(E_all(:,2),2,[]),T_ref,t_lo(:),t_hi(:));
                                        [EloTR,EhiTR] = put_error_in_time_bins_2d(reshape(E_all(:,3),2,[]),T_ref,t_lo(:),t_hi(:));
                                        [EloTL,EhiTL] = put_error_in_time_bins_2d(reshape(E_all(:,4),2,[]),T_ref,t_lo(:),t_hi(:));
                                        % place error vectors in order                                     % error in 
                                        e_verts = [EloBL EhiBL EloBR EhiBR EloTR EhiTR EloTL EhiTL];
                                        if debug_flag
                                            EloBL_T=interp1((t_lo(:)+t_hi(:))/2,EloBL,T_ref)';
                                            BLlo_Z = EloBL_T+Z_ref;
                                            
                                            plot(BLlo_Z(1,:),BLlo_Z(2,:));
                                            EhiBL_T=interp1((t_lo(:)+t_hi(:))/2,EhiBL,T_ref)';
                                            BLhi_Z = EhiBL_T+Z_ref;
                                            plot(BLhi_Z(1,:),BLhi_Z(2,:));
                                        end
                                        error_verts=[error_verts e_verts];
                                    end
                                end
                            end
                        end
                    end

                    %% summarize error
                    for i = 1:N_t_bins
                        %                             hold on;
                        verts=reshape(error_verts(i,:)',2,[]);
                        % find bounding box of the error vectors
                        [c,center_bb,h,gen] = minBoundingBox(verts);
                        if plot_flag
                            figure(2);clf;hold on;axis equal;scatter(verts(1,:),verts(2,:));
                            plot(c(1,:),c(2,:))
                        end
                        error_data_cur(1:2,i)= center_bb;
                        error_data_cur(3,i)  = h;
                        error_data_cur(4:5,i)= gen;
                    end%     toc
                    
                    error_table{v_idx,y_idx,h_idx,del_idx,vd_idx} = error_data_cur;
                    res {v_idx,y_idx,h_idx,del_idx,vd_idx,1} =Rcont;
                    res {v_idx,y_idx,h_idx,del_idx,vd_idx,2} =Rcont2;
                    res {v_idx,y_idx,h_idx,del_idx,vd_idx,3} =Rcont3;
                    if plot_flag
                        err_zono = {}; err_zono_clean = {};
                        rot_vert = rotmat(0.2)*A.footprint_vertices;
                        ft_print_gen = (max(rot_vert,[],2)-min(rot_vert,[],2))/2;
                        
                        %% viz error
                        for i = 1:N_t_bins%consider each time bin seperately
                            center_bb = error_data_cur(1:2,i);
                            h =  error_data_cur(3,i);
                            gen = error_data_cur(4:5,i);
                            len = gen(1);
                            width = gen(2);
                            ego_gen = [[cos(h)*len; sin(h)*len], [sin(-h)*width; cos(-h)*width]];
                            err_zono{i}  = zonotope([[center_bb(1);center_bb(2)], ego_gen]);
                            
                        end
                        my_rand_h = khc -h_ini_range + 2*rand*h_ini_range;
                        my_rand_vi = kvc -kvig + 2*rand*kvig;
                        my_rand_del = delta_c-del_ini_range + 2*rand*del_ini_range;
                        
                        my_rand_vd = kvd -kvdg + 2*rand*kvdg;
                        if (my_rand_vi-my_rand_vd)> 3.9
                            my_rand_vd = my_rand_vi-3.9;
                        end
                        
                        my_rand_yd = kyc -kyg + 2*rand*kyg;
                        ini_cond = [0;0;my_rand_h;my_rand_vi;my_rand_del];
                        
                        
                        
                        A.reset(ini_cond);
                        [T,U,Z] = parameterized_traj_1d_with_brake(my_rand_vd, my_rand_vi,my_rand_h,my_rand_yd);
                        A.move(t_total,T,U,Z);
                        close all;hold on;axis equal;
                        
                        for i = 1:length(Rcont)
                            %  [0;0;kvdg;kvig;h_ini_range;kyg;0]
                            z_with_err = err_zono{i}+project(zonotope_slice(Rcont{i}{1}, [3;4;5;6], [my_rand_vd;my_rand_vi;my_rand_h;my_rand_yd]), [1, 2]);
                            p3= plot(z_with_err,[1,2],'r');

                        end
                        peak_count = length(Rcont);
                        for i = 1:length(Rcont2)
                            %  [0;0;kvdg;kvig;h_ini_range;kyg;0]
                            z_with_err = err_zono{i+peak_count}+project(zonotope_slice(Rcont2{i}{1}, [3;4;5;6], [my_rand_vd;my_rand_vi;my_rand_h;my_rand_yd]), [1, 2]);
                            p3= plot(z_with_err,[1,2],'r');

                        end
                        two_sec_count = peak_count + length(Rcont2);
                        for i = 1:length(Rcont3)
                            %  [0;0;kvdg;kvig;h_ini_range;kyg;0]
                            z_with_err = err_zono{i+two_sec_count}+project(zonotope_slice(Rcont3{i}{1}, [3;4;5;6], [my_rand_vd;my_rand_vi;my_rand_h;my_rand_yd]), [1, 2]);
                            p3= plot(z_with_err,[1,2],'r');
                      
                        end
                        
                        A.plot();
                        
                        
                        %                                                          pause(3)
                        %  text(0,14, "Zonotope Uncertainty",'Color','green','FontSize',13)
                        %                         text(20,14, "Random Simulation",'Color','red','FontSize',13)
                        text(0,10,["initial heading = "+ num2str(rad2deg(khc)) + " +-0.5 deg", "initial velocity = "+num2str(kvc) + "+- 1 m/s",...
                            "desired velocity = "+ num2str(kvd)+ "+- 1 m/s","desired y = "+ num2str(kyc)+ "+- 0.5 m"])
                        %                         text(20,10,["initial heading = "+ num2str(rad2deg(my_rand_h),2) + "deg", "initial velocity = "+num2str(my_rand_vi,3) + "m/s",...
                        %                             "initial delta = "+num2str(rad2deg(my_rand_del),2)+"deg","desired velocity = "+ num2str(my_rand_vd,3)+ "m/s","desired y = "+ num2str(my_rand_yd,2)+ "m"])
                        yline(-2,'--','LineWidth',3);
                        yline(2,'--','LineWidth',3);
                        
                        axis([-5 35 -5 15]);drawnow

                    end
                    

                end
            end
        end
    end
end


%%
time = toc
% save(['zono_full_9.3_spd8.mat'],'res','error_table','v_range','kvig','v_des_range', 'kvdg','y_range','kyg','h_range','h_ini_range','del_range','del_ini_range','t_peak','dt');

%%
if debug_flag
    vi = 4 ; vf = 5;
    v0 = vi;
    v2 = vf;
    v6 = 0;
    v4 = vf/2;
    t = [0 2 4 6];
    p=polyfit(t,[v0 v2 v4 v6],4);
    t = linspace(0,6);
    ploty = polyval(p,t);
    figure()
    plot(t, ploty)
    
    y_des = -1;
    y0 = 0;
    y2 = y_des/ 2;
    y4 = y_des;
    y5 = y_des;
    y6 = y_des;
    t = [0 2 4  5 6];
    p=polyfit(t,[y0 y2 y4 y5 y6],6);
    t = linspace(0,6);
    ploty = polyval(p,t);
    figure()
    plot(t, ploty)
end
