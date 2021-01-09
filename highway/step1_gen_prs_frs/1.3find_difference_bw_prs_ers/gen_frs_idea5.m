clear; clc;close all
plot_flag = 0;
w = warning ('off','all');

v_range =2:2:4;% set the centers of the velocity, y_desired, heading and delta
%v = 1---3 , 3----5
v_des_range = v_range;
y_range = [0] ;
h_range = -0.3:0.1:0.3;
del_range = [-0.2 0  0.2 ];

%set the size of the generators
del_ini_range = 0.1; 
h_ini_range = 0.05;
kyg = 1;
kvdg = 1; % can be different for different speeds, keep it same for now
kvig = 1;

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
debug_flag = 0;


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

% set limits on the range of avaliable commanded velocity, v_cur - 4 to
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
                    my_counter= my_counter +1
                    
                    
                    kvc = v_range(v_idx); % set zonotope centers
                    kyc = y_range(y_idx);
                    
                    khc = h_range(h_idx);
                    kvd = v_range(vd_idx);
                    delta_c = del_range(del_idx);
                   
                    %% section 1 to peak
                                % ref| parameters   | time
                    num_prs = 7; %v y k_pk, kv,kh,ky, t   
                    %initial conitions are treated the same as target
                    %parameters, only difference is when slicing, initial
                    %conditions are directly sliced, but target ones are
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
                    % calculate error
                    
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
                    %     Ez_hi_cur = zeros(N_t_bins,1) ;
                    %                     Ez_hi_all = zeros(N_t_bins,N_centers) ;
                    %                     if include_ers_flag
                    %% sim hi fidelity to get error
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
                    %                     end
                    %error at time 1 %
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
                        %                         ha= {};
                        %                         for i = 1: t_total/dt
                        %                             ha{i} = [];
                        %                         end
                        %                         %                        ha = zeros(2, 0,t_total/dt);
                        %                         for i = 1:num_sims
                        % %                             rng('shuffle')
                        my_rand_h = khc -h_ini_range + 2*rand*h_ini_range;
                        my_rand_vi = kvc -kvig + 2*rand*kvig;
                        my_rand_del = delta_c-del_ini_range + 2*rand*del_ini_range;
                        
                        my_rand_vd = kvd -kvdg + 2*rand*kvdg;
                        if (my_rand_vi-my_rand_vd)> 3.9
                            my_rand_vd = my_rand_vi-3.9;
                        end
                        
                        my_rand_yd = kyc -kyg + 2*rand*kyg;
                        ini_cond = [0;0;my_rand_h;my_rand_vi;my_rand_del];
                        
                        %
                        %                         my_rand_vi = 9.776
                        %
                        %                         my_rand_vd = 9.4885
                        
                        
                        A.reset(ini_cond);
                        [T,U,Z] = parameterized_traj_1d_with_brake(my_rand_vd, my_rand_vi,my_rand_h,my_rand_yd);
                        %                             Z = [0 0 0 0 0;0 0 0 0 0;0 0 0 0 0];
                        A.move(t_total,T,U,Z);
                        close all;hold on;axis equal;
                        
                        for i = 1:length(Rcont)
                            %  [0;0;kvdg;kvig;h_ini_range;kyg;0]
                            z_with_err = err_zono{i}+project(zonotope_slice(Rcont{i}{1}, [3;4;5;6], [my_rand_vd;my_rand_vi;my_rand_h;my_rand_yd]), [1, 2]);
                            p3= plot(z_with_err,[1,2],'r');
                            %                             p3.FaceAlpha = 0.05;
                            %                         plot_arr{i}= p3;
                            %                         z_with_clean_err = err_zono_clean{i}+project(zonotope_slice(R{i}{1}, [3;4;5;6], [my_rand_vd;my_rand_vi;my_rand_h;my_rand_yd]), [1, 2])
                            %                         p4= plotFilled(z_with_clean_err,[1,2],'r');
                            %                         p4.FaceAlpha = 0.05;
                        end
                        peak_count = length(Rcont);
                        for i = 1:length(Rcont2)
                            %  [0;0;kvdg;kvig;h_ini_range;kyg;0]
                            z_with_err = err_zono{i+peak_count}+project(zonotope_slice(Rcont2{i}{1}, [3;4;5;6], [my_rand_vd;my_rand_vi;my_rand_h;my_rand_yd]), [1, 2]);
                            p3= plot(z_with_err,[1,2],'r');
                            %                             p3.FaceAlpha = 0.05;
                            %                         plot_arr{i}= p3;
                            %                         z_with_clean_err = err_zono_clean{i+peak_count}+project(zonotope_slice(Rcont2{i}{1}, [3;4;5;6], [my_rand_vd;my_rand_vi;my_rand_h;my_rand_yd]), [1, 2])
                            %                         p4= plotFilled(z_with_clean_err,[1,2],'r');
                            %                         p4.FaceAlpha = 0.05;
                        end
                        two_sec_count = peak_count + length(Rcont2);
                        for i = 1:length(Rcont3)
                            %  [0;0;kvdg;kvig;h_ini_range;kyg;0]
                            z_with_err = err_zono{i+two_sec_count}+project(zonotope_slice(Rcont3{i}{1}, [3;4;5;6], [my_rand_vd;my_rand_vi;my_rand_h;my_rand_yd]), [1, 2]);
                            p3= plot(z_with_err,[1,2],'r');
                            %                             p3.FaceAlpha = 0.05;
                            %                         plot_arr{i}= p3;
                            %                         z_with_clean_err = err_zono_clean{i+two_sec_count}+project(zonotope_slice(Rcont3{i}{1}, [3;4;5;6], [my_rand_vd;my_rand_vi;my_rand_h;my_rand_yd]), [1, 2])
                            %                         p4= plotFilled(z_with_clean_err,[1,2],'r');
                            %                         p4.FaceAlpha = 0.05;
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
                        % %                         for i = 1:length(Rcont)
                        % %                             %  [0;0;kvdg;kvig;h_ini_range;kyg;0]
                        % %                             delete(plot_arr{i});
                        % %                         end
                    end
                    
                    
                    %%
                    
                    
                    %                     res{v_idx,y_idx,h_idx,del_idx,vd_idx} = Rcont;
                    %
                    %                 for rcont_idx = 1:length(Rcont)
                    %                     res_simple{v_idx,y_idx,h_idx,vd_idx}{rcont_idx} = polygon(project(Rcont{rcont_idx}{1}+footprint_buffer,[1 2]));
                    %
                    %                 end
                end
            end
        end
    end
end
%                             figure(1);
%save all 100 sims;
%                             rng(0)

%                             for t_idx = 1:t_total/dt
%                                 t = (t_idx -1) * dt;
%                                 whichones = A.time >= t & A.time < t+dt;
%                                 %                              scatter(A.state(1,whichones),A.state(2,whichones),'MarkerEdgeColor',[rand rand rand]);
%                                 x_array = (A.state(1,whichones)) ;
%                                 h_array = A.state(3,whichones);
%                                 y_array = A.state(2,whichones);
%                                 for tim_div_idx = 1: length(x_array)
%                                     buff_vert = rotmat(h_array(tim_div_idx))*A.footprint_vertices+diag([x_array(tim_div_idx);y_array(tim_div_idx)])*ones(size(A.footprint_vertices));
%                                     ha{t_idx}= [ha{t_idx} buff_vert];
%                                 end
%                                 %         t = pi*rand(1);
%                                 %     X = [cos(t) -sin(t) ; sin(t) cos(t)]*[7 0; 0 2]*rand(2,n);
%                                 %     X = [X  20*(rand(2,1)-0.5)];  % add an outlier
%
%                                 %     tic
%
%
%
%                             end
%                        res{v_idx,y_idx,h_idx,del_idx,vd_idx} = Rcont;
%                         h=plot(A.state(1,:),A.state(2,:));
%                         ha{i}= h;
%                          A.animate();

%                     T= [0 t_brake] ;
%                     U =repmat( [0;my_rand_yd],1, length(T));
%                     Z = [0 0 0 0 0;0 0 0 0 0];
%                     A.move(t_brake,T,U,Z) ;
%                         figure(2);plot(A.time(1:end-1),diff(A.state(4,:))./diff(A.time));xline(0.5);title('acc')
%                         figure(3);clf; hold on;plot(A.time,A.state(4,:));xline(0.5);yline(0);title('spd');text(1.5, 2, "t brk ="+num2str(t_brake));
%

%                         end

%%
%                         zono_series= {};
%                         for t_idx = 1:t_total/dt
%                             [c,center,h,gen] = minBoundingBox(ha{t_idx});
%                             %     toc
%                             len = gen(1);
%                             width = gen(2);
%
%                             obs_gen = [[cos(h)*len; sin(h)*len], [sin(-h)*width; cos(-h)*width]];
%                             obs_zono = zonotope([[center(1);center(2)], obs_gen]);
%                             zono_series{t_idx} = obs_zono;
%                             if plot_flag
%
%                             hpts  = plot(ha{t_idx}(1,:),ha{t_idx}(2,:),'.');
%                             hcont = plot(c(1,[1:end 1]),c(2,[1:end 1]),'r');
%                              p_obs = plotFilled(obs_zono, [1, 2], 'r');
%                             p_obs.FaceAlpha = 0.3;
%                             p_obs.EdgeAlpha = 0.3;
%                             delete(hpts);delete(hcont);delete(p_obs);
%                             %%
%                             end
% %                             pause(0.5)
%
%
%
%
%                         end
%                         res{v_idx,y_idx,h_idx, del_idx ,vd_idx } = zono_series;
%                         for i = 1:100
%                         delete(ha{i})
%                         end
%
%  figure(2); clf;hold on; title('pos error propogation')

%                         for i=1:length(simRes.t)
%                             plot(simRes.x{i}(:,1),simRes.x{i}(:,2),'b');
%                         end

%%
%
% %                          for i = 1:length(Rcont_stop)
% %                             Rcont_stop{i}{1}=Rcont_stop{i}{1}+footprint_buffer;
% %                             p2 = plotFilled(Rcont_stop{i}{1},[1, 2], 'r');
% %                             p2.FaceAlpha = 0.2;
% %                             %                                     if i == 1
% %                             %                                         p = plotFilled(Rcont{i}{1}, [1, 2], 'r');
% %                             %                                         p.FaceAlpha = 0.3;
% %                             %
% %                             %                                     elseif i == length(Rcont)
% %                             %                                         p = plotFilled(Rcont{i}{1}, [1, 2], 'b');
% %                             %                                         p.FaceAlpha = 0.3;
% %                             %
% %                             %                                     else
% %                             %                                         p = plot(Rcont{i}{1}, [1, 2], 'g');
% %                             %
% %                             %                                     end
% %                          end
%                         %                                                     pause(2)
%
%                         for i = 1:length(Rcont)
%                             p3=plotFilled(zonotope_slice(Rcont{i}{1}, [6;7], [my_rand_vd;my_rand_yd]), [1, 2], 'b');
%                             p3.FaceAlpha = 0.2;
%                         end
%                         v_array = [];
%                         for i = 1:length(Rcont)
%                             p3=plotFilled(zonotope_slice(Rcont{i}{1}, [6;7;8;9;10], [my_rand_vd;my_rand_yd;my_rand_h;my_rand_vi;my_rand_del]), [1, 2], 'r');
%                             p3.FaceAlpha = 0.2;
%                             p4=deleteAligned(project(zonotope_slice(Rcont{i}{1}, [6;7;8;9;10], [my_rand_vd;my_rand_yd;my_rand_h;my_rand_vi;my_rand_del]), [4]));
%                             v_array = [v_array; my_rand_vi+p4.center];
%                             if(my_rand_vi +p4.center)< 0.1
%                                  break
%                             end
%                         end
%
% %                         for i = 1:length(Rcont_stop)
% %                             z_slice = zonotope_slice(Rcont_stop{i}{1}, [6;7;8;9;10], [my_rand_vd;my_rand_yd;my_rand_h;my_rand_vi;my_rand_del]);
% %                             z_v=deleteAligned(project(z_slice,[4]));
% %                             if(my_rand_vi +z_v.center)< 0
% %                                 break
% %                             end
% %                             p3=plotFilled(z_slice, [1, 2], 'y');
% %                             p3.FaceAlpha = 0.1;
% %
% %                         end
%
%                         % %                                                     pause(2)
% %                         A.animate()
%                         figure(4);
%                         plot(dt:dt:t_peak,v_array,dt:dt:t_peak-dt,diff(v_array)./diff(dt:dt:t_peak))
%                         figure(2);plot(A.time(1:end-1),diff(A.state(4,:))./diff(A.time));xline(0.5);title('acc')
%                         figure(3);clf; hold on;plot(A.time,A.state(4,:));xline(0.5);yline(0);title('spd');text(1.5, 2, "t brk ="+num2str(t_brake));
%                         %plot(A.time,(A.state(4,:) - 4) .* (-(A.time-0.5)+5)./5);
% %                         vq1 = interp1(T,U(1,:),A.time); plot(A.time, vq1);
%                         figure(1)
%
%                     kvc = 8.5; kvd = 9.5;

% 22 + 2*(1- 2);
%                 if kvd < 22 || kvd > 32
% %                     disp(['kvd exceeed bound, skipping kvd = ', num2str(kvd)]);
%                     continue;
%                 end
%                     options = struct;
%                     options.tStart = 0;
%                     options.tFinal = t_peak;
%                     options.timeStep=dt; %time step size for reachable set computation
%                     options.taylorTerms=3; %number of taylor terms for reachable sets
%                     options.zonotopeOrder= 3; %zonotope order... increase this for more complicated systems.
%                     options.maxError =9e100*ones(num_zono_state, 1); % this controls splitting, set it high to avoid splitting
%                     options.verbose = 1;
%
%                     options.uTrans = 0; % set input = 0
%                     options.U = zonotope([0,0 ]);
%                     options.advancedLinErrorComp = 0;
%                     options.tensorOrder = 1;
%                     options.reductionInterval = inf;
%                     options.reductionTechnique = 'girard';
%                     %%
%                     sys = nonlinearSys(num_zono_state, 1, @closed_loop_highwaySS_together, options);
%                     % y ref has discrete values also
%
%                     % if kvc-v_ini_range < 22
%                     %     v_ini_tol = kvc;
%                     % elseif kvc+v_ini_range>2
%                     %     v_ini_range = 2-kvc;
%                     % else
%                     %     v_ini_tol = v_ini_range;
%                     % end
%
%                     %            X y delh delvx  deldelta; vx_des y_des; h0 vx0 delta0;time
%                     options.x0 = [0;0;0;0;0; kvd;kyc; khc;kvc; delta_c;0];
%
%                     options.R0 = zonotope([options.x0,diag([0;0;0;0;0;kvdg;kyg;h_ini_range;kvig;del_ini_range;0])]);
%
%
%

%                     Rcont = reach(sys, options);
%                     sys = nonlinearSys(num_zono_state, 1, @closed_loop_highwaySS, options);
%                     runs = 20;
%                     fracV = 0.5;
%                     fracI = 0.5;
%                     changes = 6;
%                     simRes = simulate_random(sys, options, runs, fracV, fracI, changes);
%
%                 %% computing failsafe trajectory
%                      p4=deleteAligned(project(Rcont{end}{1},[4]));
%                                   %cur_max_spd = ini max spd + min change.
%                      t_brake =round(( kvc+kvg - (p4.center - p4.generators))/2.6,1);
%                      if t_brake < 0.1
%                         t_brake = 0.1;
%                      end
%
% %                     if kvc > 4.5
%                         options.timeStep = dt_stop;
%                         options.tFinal = t_brake;
%                         options.R0 = Rcont{end}{1};
%                         options.x0 = center(options.R0)- [0; 0; 0; 0; 0;0;0;0;0;0;t_peak];
%                         options.U = zonotope([t_brake, 0]);
%                         sys = nonlinearSys(num_zono_state, 1, @closed_loop_highwaySS_brake, options);
%                         Rcont_stop = reach(sys, options);
% % %                     else
%                         options.tFinal = t_brake;
%                         options.R0 = Rcont{end}{1};
%                         options.x0 = center(options.R0)- [0; 0; 0; 0; 0;0;0;0;0;0;t_peak];
%                         sys = nonlinearSys(num_zono_state, 1, @closed_loop_highwaySS_brake, options);
%                         Rcont_stop = reach(sys, options);
%                     end

% % % %                 tComp = toc;
% % % %                 disp(['computation time of reachable set: ', num2str(tComp)]);
% % % %                 %%
% % % %                 slicek1 = 0.8
% % % %                 slicek2 = 0.2
% % % %                 slice_dim = [1;2]; % rows corresponding to parameters
% % % %                 slice_parameter = [slicek1; slicek2]; %dx dy
% % % %                 % [t,y]=ode45(@(t,x)closed_loop_turtlebot(t,x,[0]),[0  t_total],[options.x0(1:4);slice_parameter]);
% % % %
% % % %                 %%



%                 res_simple{v_idx,w_idx,h_idx,vd_idx} = cell(t_total/dt,1);

%%
time = toc
save(['zono_full_9.3_spd8.mat'],'res','error_table','v_range','kvig','v_des_range', 'kvdg','y_range','kyg','h_range','h_ini_range','del_range','del_ini_range','t_peak','dt');

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
% save(['Zonotope_simple_2sec_6.4.mat'],'res_simple','-v7.3');
% for i = 1:length(Rcont)
%     p2 = plotFilled(zonotope_slice(Rcont{i}{1}, slice_dim, slice_parameter), [1, 2], 'b');
%     p2.FaceAlpha = 0.25;
% end
%%
% figure(2);clf;
%  xlabel("heading");ylabel('velocity');xlim([-1 1]);ylim([-1 1])
% figure(3)
%  xlabel('x_error');ylabel ('y_error');xlim([-1 1]);ylim([-1 1])