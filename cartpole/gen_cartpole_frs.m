clear all;
plot_flag = 0;
% A= cartpole_agent;
num_range = 5;
N_theta =4;% initial
N_theta_dot = 3;% initial
N_x_dot = num_range;% initial & parameter
N_ka    = num_range;% initial & parameter
N_kpeak = num_range;% selectable &  parameter

theta_start = -pi; theta_end = pi/2;
theta_dot_start= -2; theta_dot_end = 2; % no idea what the range of this should be :(
v_start = -5; v_end = 5;
ka_start = -15; ka_end = 15;
kpeak_start = -5; kpeak_end = 5;


% sample the parameter space
theta_a = linspace(theta_start,theta_end,N_theta);
theta_dot_a = linspace(theta_dot_start,theta_dot_end,N_theta_dot);
vi = linspace(v_start,v_end,11);
ka = linspace(ka_start,ka_end,N_ka);
kpeak = vi;%/linspace(kpeak_start,kpeak_end,N_kpeak);

kvg = abs(vi(1)-vi(2))/2;
kag = abs(ka(1)-ka(2))/2;
kpkg = abs(kpeak(1)-kpeak(2))/2;
Tg = abs(theta_a(1)-theta_a(2))/2;
dTg = abs(theta_dot_a(1)-theta_dot_a(2))/2;

Tf = 0.3;
Tpk = 0.1;
dt = 0.01

t_lo= 0:dt:Tf-dt; % time sections to contain error data
t_hi= dt:dt:Tf;
N_t_bins = length(t_lo);

A = cartpole_agent;

res=cell(length(vi),length(ka),3,2); % each section has at most 3 possible desired velocity
error_table = cell(length(vi),length(ka),3,length(theta_a),length(theta_dot_a));

if plot_flag
    figure('Renderer', 'painters', 'Position', [10 10 800*1.5 400*1.5])
end
tic
for v_idx = 1:length(vi)
    v_idx
    for vd_idx = 1:3
        vd_real_idx = v_idx+vd_idx-2;
        if  vd_real_idx < 1 || vd_real_idx > length(kpeak)
            continue;
        end
        for ka_idx = 1:length(ka)
            
            kvc = vi(v_idx);%inital velocity param
            kpkc = kpeak(vd_real_idx);
            kac = ka(ka_idx);

            num_prs = 5; %v; k_pk; kv;ka; t1  % the first element v is not the actual v but delta v from kv
                    %
            options = struct;
            options.tStart = 0;
            options.tFinal = Tpk;
            options.timeStep=dt; %time step size for reachable set computation
            options.taylorTerms=3; %number of taylor terms for reachable sets
            options.zonotopeOrder= 3; %zonotope order... increase this for more complicated systems.
            options.maxError =9e100*ones(num_prs, 1); % this controls splitting, set it high to avoid splitting
            options.verbose = 0;
                    
            options.uTrans = 0; % set input = 0
            options.U = zonotope([0,0]);
            options.advancedLinErrorComp = 0;
            options.tensorOrder = 1;
            options.reductionInterval = inf;
            options.reductionTechnique = 'girard';
                    
            sys = nonlinearSys(num_prs, 1, @cartpole_toPeak, options);
            options.x0 = [0;kpkc;kvc;kac;0];
                                 %v y v_desired@2sec, vini, hini, y_desired@4sec, t
            options.R0 = zonotope([options.x0,diag([0;kpkg;kvg;kag;0])]);
            Rcont = reach(sys, options);
                    
            %% section 2 stop 
            % set final condition of previous zonotope as initial 
            % condition of this one, also decrease time
            options.R0 = Rcont{end}{1} - [0;0;0;0;Tpk]; %initial state for reachability analysis
            options.x0 = center(options.R0);
            options.tFinal = Tf -Tpk;
            sys = nonlinearSys(num_prs, 1, @cartpole_toStop, options);
            Rcont2 = reach(sys, options);
            
            res{v_idx,ka_idx,vd_idx,1} = Rcont;
            res{v_idx,ka_idx,vd_idx,2} = Rcont2;
            
            zono_one_all = [Rcont;Rcont2];
            n = length(zono_one_all);
%             %calculate error
            for theta_idx = 1:length(theta_a)
                for theta_dot_idx = 1:length(theta_dot_a)
                    Tc = theta_a(theta_idx);
                    dTc = theta_dot_a(theta_dot_idx);
                    
                    
                    Ex_lo = inf*zeros(N_t_bins,1) ;
                    Ex_hi = -inf*zeros(N_t_bins,1) ;
                    
                    num_sims = 3;
                    for kv_sim = linspace(kvc-kvg, kvc+kvg, num_sims)
                        for kpk_sim = linspace(kpkc-kpkg,kpkc+kpkg,num_sims)
                            for ka_sim = linspace(kac-kag, kac+kag, num_sims)
                                [T_ref,U_ref,Z_ref] =parameterized_cartpole_traj(kpk_sim,kv_sim,ka_sim);
                                if plot_flag
                                figure(1); clf; hold on;
                                sh_prs = scatter([],[],30);
                                sh_prss = scatter([],[],10);
                                plot(T_ref, U_ref, T_ref, Z_ref);
                                for t_idx = 1:n
                                    zono_one = zono_one_all{t_idx}{1};
                                    zono_one_sliced = zonotope_slice(zono_one, [2;3;4], [kpk_sim;kv_sim;ka_sim]);
                                    zono_prss_vert = deleteAligned(project(zono_one_sliced, 1));
                                    zono_prs_vert = deleteAligned(project(zono_one, 1));
                                    pshi = center(zono_prss_vert) +generators(zono_prss_vert);
                                    pslo = center(zono_prss_vert) -generators(zono_prss_vert);
                                    t_now = 0.3/n*t_idx;
                                    sh_prss.XData = [sh_prss.XData t_now t_now];
                                    sh_prss.YData = [sh_prss.YData pslo pshi];
                                    
                                    phi = center(zono_prs_vert) +generators(zono_prs_vert);
                                    plo = center(zono_prs_vert) -generators(zono_prs_vert);
                                    sh_prs.XData = [sh_prs.XData t_now t_now];
                                    sh_prs.YData = [sh_prs.YData plo phi];
                                end
                                end
%                                     prs_zono = zonotope_slice(zono_one, 2, K);
                                for Tsim =linspace(Tc-Tg, Tc+Tg, num_sims)
                                    for dTsim = linspace(dTc-dTg, dTc+dTg, num_sims)
                                               % x  dx     T    dT ka
%                                          A.reset([0;kv_sim;Tsim;dTsim;ka_sim]);
                                        [t, y] = ode45(@(t,y)A.cartpole_dynamic_new(t,y,T_ref,U_ref,Z_ref),[0 Tf],[0;kv_sim;Tsim;dTsim]);
                                        if plot_flag
                                            plot(T_ref,U_ref, T_ref, Z_ref,t,y(:,1),t,y(:,2));
                                            legend('vnom','xnom','x','v')
                                            drawnow
                                            pause(0.1)
                                        end
                                        Z = match_trajectories(T_ref,t,y');
                                        [Elo,Ehi]=put_error_in_time_bins_1d(Z-Z_ref,T_ref,t_lo(:),t_hi(:));
                                        Ex_lo = min(Elo, Ex_lo);
                                        Ex_hi = min(Ehi, Ex_hi);
                                        
                                    end
                                end
                            end
                        end
                    end
                    
                   
                    
                    %place planning reach set inside data structure
                    %simulate range in the zonotope to genreate ERS
                    error_table{v_idx,ka_idx,vd_idx,theta_idx,theta_dot_idx}=[Ex_lo Ex_hi];
                    
                end
            end
            
            

            
        end
    end
end
%%
toc
save('cartpole_frs_9_22.mat','res','error_table','theta_a','theta_dot_a','vi','ka','kpeak','kvg','kag','kpkg','Tg','dTg','Tf','Tpk','dt');
% vd_idx_max_a = zeros(length(v_range),1);
% vd_idx_min_a = zeros(length(v_range),1);
