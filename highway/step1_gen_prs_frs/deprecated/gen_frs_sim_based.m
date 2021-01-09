clear; clc;close all
plot_flag = 0;

num_zono_state = 11;
v_range = 0.5:1:9.5;%[9.5 8.5 7.5 6.5 5.5 4.5 3.5 2.5 1.5 0.5];
num_v_des =length(v_range);
% y_range_first_half = 0.2:0.4:1.8;
% y_range_second_half = -1.8:0.4:0;
y_range =  [-2 2] ; % since may be in the middle of a lane change behavior, actually
h_range = [-0.15 -0.05 0.05 0.15];
% y_range = [2 6 10];
% y_range = [y_range_second_half y_range_first_half];
del_range = [-0.05 0.05];
del_ini_range = 0.05; %0.07; %from  experiment
h_ini_range = 0.05;
kyg = 2;
kvdg = 0.5;
kvig = 0.5;
num_sims =100;

rot_max = 0.2;
rot_len = cos(rot_max)* 4.8 + sin(rot_max)*2;
rot_width = sin(rot_max)*4.8 + cos(rot_max)*2;


footprint_buffer = zonotope([[0;0;0;0;0;0;0;0;0;0;0],diag([rot_len/2;rot_width/2;0;0;0;0;0;0;0;0;0 ])]);

v_des_range = v_range;
t_peak = 0.5;
dt = 0.2;
dt_stop = 0.1;
yr = 1:length(y_range);
vr = 1:length(v_range);
hr = 1:length(h_range);
deltar = 1:length(del_range);

%t_brake =4;
t_total = 5;
t_brake= t_total - t_peak;
tic
%%
res=cell(length(v_range),length(y_range),length(h_range),length(del_range),num_v_des);
res_stop=cell(length(v_range),length(y_range),length(h_range),length(del_range),num_v_des);
% res_simple=cell(length(v_range),length(y_range),length(h_range),length(del_range),num_v_des);
A = highway_cruising_agent;
if plot_flag
    
    figure('Renderer', 'painters', 'Position', [10 10 800*1.5 400*1.5])
end

vd_idx_max_a = zeros(length(v_range),1);
vd_idx_min_a = zeros(length(v_range),1);
pts = [];
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
figure(3);
scatter(pts(1,:),pts(2,:));
%
% for h_idx = 1: hr
my_counter = 0;
rng('shuffle')
for v_idx = vr
    for h_idx =hr
        for vd_idx =vd_idx_min_a(v_idx) :vd_idx_max_a(v_idx)
            for del_idx = deltar
                for y_idx = yr
                    my_counter= my_counter +1
                    
                    kvc = v_range(v_idx);
                    kyc = y_range(y_idx);
                    
                    khc = h_range(h_idx);
                    kvd = v_range(vd_idx);
                    delta_c = del_range(del_idx);
%                     figure(1);hold on
                        %                         clf;hold on;
                        % for i = 1:length(Rcont)
                        %     p2 = plotFilled(zonotope_slice(Rcont{i}{1}, slice_dim, slice_parameter), [3, 4], 'b');
                        %     p2.FaceAlpha = 0.25;
                        % end
                        % plot(y(:,3),y(:,4),'r');
                        
                      
                        ha= {};
                        for i = 1: t_total/dt
                            ha{i} = [];
                        end
                        %                        ha = zeros(2, 0,t_total/dt);
                        for i = 1:num_sims
%                             rng('shuffle')
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
                            T= [0         t_peak            t_total] ;
                            U=[my_rand_vd my_rand_vd                      0 ;
                                my_rand_yd my_rand_yd      my_rand_yd ];
                            
                            Z = [0 0 0 0 0;0 0 0 0 0;0 0 0 0 0];
                            A.move(t_total,T,U,Z) ;
                            
%                             figure(1);
                            %save all 100 sims;
%                             rng(0)
                            
                            for t_idx = 1:t_total/dt
                                t = (t_idx -1) * dt;
                                whichones = A.time >= t & A.time < t+dt;
                                %                              scatter(A.state(1,whichones),A.state(2,whichones),'MarkerEdgeColor',[rand rand rand]);
                                x_array = (A.state(1,whichones)) ;
                                h_array = A.state(3,whichones);
                                y_array = A.state(2,whichones);
                                for tim_div_idx = 1: length(x_array)
                                    buff_vert = rotmat(h_array(tim_div_idx))*A.footprint_vertices+diag([x_array(tim_div_idx);y_array(tim_div_idx)])*ones(size(A.footprint_vertices));
                                    ha{t_idx}= [ha{t_idx} buff_vert];
                                end
                                %         t = pi*rand(1);
                                %     X = [cos(t) -sin(t) ; sin(t) cos(t)]*[7 0; 0 2]*rand(2,n);
                                %     X = [X  20*(rand(2,1)-0.5)];  % add an outlier
                                
                                %     tic
                                
                                
                                
                            end
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
                            
                        end
                        
                        %%
                        zono_series= {};
                        for t_idx = 1:t_total/dt
                            [c,center,h,gen] = minBoundingBox(ha{t_idx});
                            %     toc
                            len = gen(1);
                            width = gen(2);
                            
                            obs_gen = [[cos(h)*len; sin(h)*len], [sin(-h)*width; cos(-h)*width]];
                            obs_zono = zonotope([[center(1);center(2)], obs_gen]);
                            zono_series{t_idx} = obs_zono;
                            if plot_flag
                                
                            hpts  = plot(ha{t_idx}(1,:),ha{t_idx}(2,:),'.');
                            hcont = plot(c(1,[1:end 1]),c(2,[1:end 1]),'r');
                             p_obs = plotFilled(obs_zono, [1, 2], 'r');
                            p_obs.FaceAlpha = 0.3;
                            p_obs.EdgeAlpha = 0.3;
                            delete(hpts);delete(hcont);delete(p_obs);
                            %%
                            end
%                             pause(0.5)
                            
                            
                            
                            
                        end
                        res{v_idx,y_idx,h_idx, del_idx ,vd_idx } = zono_series;
                        %                         for i = 1:100
                        %                         delete(ha{i})
                        %                         end
                        %
                        %  figure(2); clf;hold on; title('pos error propogation')
                        
                        %                         for i=1:length(simRes.t)
                        %                             plot(simRes.x{i}(:,1),simRes.x{i}(:,2),'b');
                        %                         end
                        
                        %                         for i = 1:length(Rcont)
                        %                             figure(1)
                        %                             Rcont{i}{1}=Rcont{i}{1}+footprint_buffer;
                        %                             p2 = plotFilled(Rcont{i}{1},[1, 2], 'g');
                        %                             p2.FaceAlpha = 0.2;
                        %
                        %                             p4 = project(Rcont{i}{1},[4]);
                        %                             p4 = deleteAligned(p4);
                        %
                        %
                        %                             %                                     if i == 1
                        %                             %                                         p = plotFilled(Rcont{i}{1}, [1, 2], 'r');
                        %                             %                                         p.FaceAlpha = 0.3;
                        %                             %
                        %                             %                                     elseif i == length(Rcont)
                        %                             %                                         p = plotFilled(Rcont{i}{1}, [1, 2], 'b');
                        %                             %                                         p.FaceAlpha = 0.3;
                        %                             %
                        %                             %                                     else
                        %                             %                                         p = plot(Rcont{i}{1}, [1, 2], 'g');
                        %                             %
                        %                             %                                     end
                        %                         end
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
                        % %                         for i = 1:length(Rcont)
                        % %                             p3=plotFilled(zonotope_slice(Rcont{i}{1}, [6;7], [my_rand_vd;my_rand_yd]), [1, 2], 'b');
                        % %                             p3.FaceAlpha = 0.2;
                        % %                         end
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
                    if plot_flag
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
                    %                     res{v_idx,y_idx,h_idx,del_idx,vd_idx} = Rcont;
                    %                     res_stop{v_idx,y_idx,h_idx,del_idx,vd_idx} =Rcont_stop;
                    %                 for rcont_idx = 1:length(Rcont)
                    %                     res_simple{v_idx,y_idx,h_idx,vd_idx}{rcont_idx} = polygon(project(Rcont{rcont_idx}{1}+footprint_buffer,[1 2]));
                    %
                    %                 end
                end
            end
        end
    end
end


%                 res_simple{v_idx,w_idx,h_idx,vd_idx} = cell(t_total/dt,1);

%%
time = toc
del_ini_range = 0.05; %0.07; %from  experiment
h_ini_range = 0.05;
kyg = 2;
kvdg = 0.5;
kvig = 0.5;
save(['zonotpe_full_6.19.mat'], 'res','v_range','kvig','v_des_range', 'kvdg','y_range','kyg','h_range','h_ini_range','del_range','del_ini_range','t_peak','dt','-v7.3');
save(['zonotpe_stop_6.12.mat'], 'res_stop','v_range', 'y_range','h_range','del_range','t_brake','dt','-v7.3');

%%
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