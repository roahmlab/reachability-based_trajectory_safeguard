run_quadrotor_testing
%similar to car, can be used to plot smooth quadrotor motion, somewhat deprecated, have not tried in a while
close all;


%%
close all;

multi_agent_flag= true;
plot_ref_flag = 1;
plot_proposed = 1;
v = VideoWriter('drone_sim');
v.FrameRate=10;
open(v);
agent_list =[]; 

% sim_summary = load('sim_summary_5-41_15-48-23.125.mat');
sim_summary = load('sim_summary_5-41_21-24-29.435.mat');%short Z
% sim_summary = load('sim_summary_0-41_21-31-27.861.mat');%hover N
% sim_summary = load('sim_summary_1-41_21-42-37.102.mat');% crash N
offset = 10;
A.state = sim_summary.agent_info.state(:,1:end-offset);
A.attitude = sim_summary.agent_info.attitude(:,:,1:end-offset);
A.time =  sim_summary.agent_info.time(:,1:end-offset);
fh= figure('Position',[1 1 1920 1080/2])
clf;axis equal;hold on;
W.N_obstacles = length(sim_summary.world_info.obstacles);
W.obstacles =sim_summary.world_info.obstacles;
clear_plot_data(W);
W.plot();
A.plot();

if multi_agent_flag
    xlim([0 25])
    ylim([-2.5 2.5])
    zlim([-2.5 2.5])
    view(-10,30)
    num_draw = A.time (end)/0.5;
    A_arr = []
    for i =1:num_draw
        A_time = quadrotor_agent;
        A_time.plot_body_opacity = i/num_draw /2;
        A_time.time = A.time(:,1:end-offset);
        A_time.state = A.state(:,1:end-offset);
        A_time.attitude = A.attitude(:,:,1:end-offset);
        if rem(i,1)==0
            A_time.plot_at_time(0.5*i);
        end
    end
end


for i = 0:0.08:(A.time(end))%length(sim_summary.t_real_start_arr)
    
    %         plot_time = sim_summary.t_real_start_arr(i)+2;
    plot_time =i;
    real_time_arr = find(plot_time>=sim_summary.t_real_start_arr);

    real_ref_idx = (real_time_arr(end)-1)*3-2;
    if plot_proposed
        prop_time_arr = find(plot_time>=sim_summary.t_proposed_start_arr);
        prop_ref_idx = (prop_time_arr(end-1)-1)*3+1;
    end
    
    %%
    if plot_time <= A.time(end)
        if multi_agent_flag
            %                 if k > (num_drawing_steps -4)
            opa = plot_time/sim_summary.t_real_start_arr(end);
            opa=opa/4 ;
            agent_list(k).plot_footprint_opacity = opa;
            agent_list(k).plot_footprint_edge_opacity = opa*2;
            agent_list(k).plot_at_time(plot_time);
            %                 end
            k = k+1;
        else
            A.plot_at_time(plot_time);
        end
        real_time_arr = find(plot_time>=sim_summary.t_real_start_arr);
%         real_ref_idx = (real_time_arr(end)-1)*3+1;
        if plot_ref_flag
            plot_ref3(AH,[],sim_summary.ref_Z(real_ref_idx:real_ref_idx+2,:));
            if ~isempty(sim_summary.proposed_ref_Z)
                plot_ref3(AH,sim_summary.ref_Z(real_ref_idx:real_ref_idx+2,:),sim_summary.proposed_ref_Z(prop_ref_idx:prop_ref_idx+2,:))
            end
            
        end
        
    end
    
    text()
    drawnow
    frame = getframe(fh) ;
    writeVideo(v,frame);
    %     im = frame2im(frame);
    %     [imind,cm] = rgb2ind(im,256);
    
    %     if start_gif
    %         imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
    %             'DelayTime',frame_rate) ;
    %         start_gif = false ;
    %     else
    %         imwrite(imind,cm,filename,'gif','WriteMode','append',...
    %             'DelayTime',frame_rate) ;
    %     end
end
v.close();

  set(gca,'FontSize',20)
h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'drone_sim.pdf','-dpdf','-r0')