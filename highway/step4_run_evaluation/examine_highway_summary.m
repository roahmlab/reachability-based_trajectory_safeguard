run_highway_testing
close all;
compare = 0;
sim_summary=load('sim_summary_1-6_14-05-01.945.mat');
A.state = sim_summary.agent_info.state;
A.time =  sim_summary.agent_info.time;
W.obstacles_seen = sim_summary.world_info.obstacles;
AH_plot_flag = 1;
multi_agent_flag = 0;
plot_ref_flag = 1;
%%
if compare
    sim_summary2=load('sim_summary_1-10_13-02-42.670.mat');
    
    W = dynamic_car_world('bounds',bounds,'buffer',world_buffer,...
        'verbose',verbose_level,'goal_radius',goal_radius) ;
    A = highway_cruising_agent;
    % LoggedSignals=LoggedSignals.LoggedSignals;
    % LoggedSignals2=load('./savedAgents_D6Z/D6Z_sim_summary/sim_summary_5-2513_07-24-56.816.mat');
    % LoggedSignals2= LoggedSignals2.LoggedSignals;
    
    % W = dynamic_car_world('bounds',bounds,'buffer',world_buffer,...
    %     'verbose',verbose_level,'goal_radius',goal_radius) ;
    % A = highway_cruising_agent;
    
    A.state = sim_summary.agent_info.state;
    A.time =  sim_summary.agent_info.time;
    W.obstacles_seen = sim_summary.world_info.obstacles;
    AH = highwayAgentHelper(A,'zono_full_9.3_spd8.mat','t_move',t_move,'t_failsafe_move',t_failsafe_move,'eps',0.001,'verbose',verbose_level,'plot_flag',1);
    
    
    
    
    W2 = dynamic_car_world('bounds',bounds,'buffer',world_buffer,...
        'verbose',verbose_level,'goal_radius',goal_radius) ;
    A2 = highway_cruising_agent;
    % LoggedSignals=LoggedSignals.LoggedSignals;
    % LoggedSignals2=load('./savedAgents_D6Z/D6Z_sim_summary/sim_summary_5-2513_07-24-56.816.mat');
    % LoggedSignals2= LoggedSignals2.LoggedSignals;
    
    % W = dynamic_car_world('bounds',bounds,'buffer',world_buffer,...
    %     'verbose',verbose_level,'goal_radius',goal_radius) ;
    % A = highway_cruising_agent;
    
    A2.state = sim_summary2.agent_info.state;
    A2.time =  sim_summary2.agent_info.time;
    W2.obstacles_seen = sim_summary2.world_info.obstacles;
    AH2 = highwayAgentHelper(A2,'zono_full_9.3_spd8.mat','t_move',t_move,'t_failsafe_move',t_failsafe_move,'eps',0.001,'verbose',verbose_level,'plot_flag',1);
    
    fh=figure('Name','SimCompare','Position',[1 1 1000 600]); clf;
    subplot (2,1,1);hold on;
    a = gca;
    
    W.plot();
    subplot (2,1,2);hold on; axis equal;
    W2.plot();
    filename = 'animation.gif';start_gif= true;frame_rate=0;
    [max_timestep,which_agent_longer]=max([length(A.time) length(A2.time)])
    
    for j = linspace(1,max_timestep,5)
        
        i = round(j);
        subplot (2,1,1);
        if which_agent_longer == 1
            plot_time = A.time(i);
        else
            plot_time = A2.time(i);
        end
        
        real_time_arr = find(plot_time>=sim_summary.t_real_start_arr);
        real_ref_idx = (real_time_arr(end)-1)*2+1;
        prop_time_arr = find(plot_time>=sim_summary.t_proposed_start_arr);
        prop_ref_idx = (prop_time_arr(end)-1)*2+1;
        
        %%
        if plot_time <= A.time(end)
            A.plot_at_time(plot_time);
            real_time_arr = find(plot_time>=sim_summary.t_real_start_arr);
            real_ref_idx = (real_time_arr(end)-1)*2+1;
            if plot_ref_flag
                plot_ref(AH,[],sim_summary.ref_Z(real_ref_idx:real_ref_idx+1,:));
                plot_ref(AH,sim_summary.ref_Z(real_ref_idx:real_ref_idx+1,:),sim_summary.proposed_ref_Z(prop_ref_idx:prop_ref_idx+1,:))
            end
        end
        %%
        subplot (2,1,2);
        time_arr2 = find(plot_time>=sim_summary2.t_real_start_arr);
        if plot_time <= A2.time(end)
            A2.plot_at_time(plot_time);
            real_time_arr = find(plot_time>=sim_summary2.t_real_start_arr);
            real_ref_idx = (real_time_arr(end)-1)*2+1;
            if plot_ref_flag
                plot_ref(AH2,[],sim_summary2.ref_Z(real_ref_idx:real_ref_idx+1,:));
            end
        end
        
        
        drawnow
        frame = getframe(fh) ;
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        if start_gif
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
                'DelayTime',frame_rate) ;
            start_gif = false ;
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append',...
                'DelayTime',frame_rate) ;
        end
    end
    
else
    %%
    height = (bounds(4)-bounds(3))/(bounds(2)-bounds(1)) *1920;
    
    fh=figure('Name','SimCompare','Position',[1 1 1920 300]); clf;
    hold on; axis equal; W.plot();
    v = VideoWriter('car_sim');
    open(v);
    filename = 'animation.gif';start_gif= true;frame_rate=0;
    A = highway_cruising_agent;
    % LoggedSignals=LoggedSignals.LoggedSignals;
    % LoggedSignals2=load('./savedAgents_D6Z/D6Z_sim_summary/sim_summary_5-2513_07-24-56.816.mat');
    % LoggedSignals2= LoggedSignals2.LoggedSignals;
    
    % W = dynamic_car_world('bounds',bounds,'buffer',world_buffer,...
    %     'verbose',verbose_level,'goal_radius',goal_radius) ;
    % A = highway_cruising_agent;
    
    A.state = sim_summary.agent_info.state;
    A.time =  sim_summary.agent_info.time;
    if multi_agent_flag
        agent_list = []
        for  i = 1:40
            A_time = highway_cruising_agent;
            A_time.state = sim_summary.agent_info.state;
            A_time.time =  sim_summary.agent_info.time;
            agent_list = [agent_list; A_time];
        end
    end
    k = 1;
    %     sim_length = sim_summary.agent_info.time(end)/2;
    
    dt = 1;
    if multi_agent_flag
        dt = 2;
    end
    z_t = match_trajectories(0,A.time,A.state) ;
    
    th = text(double(fh.Children.XLim(2)-30),2,['Velocity = ',num2str(z_t(4))],'FontSize',20);
    for i = 0:dt:(A.time(end))%length(sim_summary.t_real_start_arr)
        
        %         plot_time = sim_summary.t_real_start_arr(i)+2;
        plot_time =i;
        real_time_arr = find(plot_time>=sim_summary.t_real_start_arr);
        real_ref_idx = (real_time_arr(end)-1)*2+1;
        if plot_ref_flag
            
            prop_time_arr = find(plot_time>=sim_summary.t_proposed_start_arr);
            if ~isempty(sim_summary.proposed_ref_Z)
                prop_ref_idx = (prop_time_arr(end)-1)*2+1;
            end
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
                z_t = match_trajectories(plot_time,A.time,A.state) ;
                th.Position=[(fh.Children.XLim(2)-30),2,0];
                th.String= ['Velocity = ',num2str(z_t(4))];
            end
            real_time_arr = find(plot_time>=sim_summary.t_real_start_arr);
            real_ref_idx = (real_time_arr(end)-1)*2+1;
            if plot_ref_flag
                plot_ref(AH,[],sim_summary.ref_Z(real_ref_idx:real_ref_idx+1,:));
                if ~isempty(sim_summary.proposed_ref_Z)
                    plot_ref(AH,sim_summary.ref_Z(real_ref_idx:real_ref_idx+1,:),sim_summary.proposed_ref_Z(prop_ref_idx:prop_ref_idx+1,:))
                end
            end
            
        end
        
        text()
        drawnow
        frame = getframe(fh) ;
        writeVideo(v,frame);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        if start_gif
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
                'DelayTime',frame_rate) ;
            start_gif = false ;
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append',...
                'DelayTime',frame_rate) ;
        end
    end
    close(v);
    xlim ([bounds(1),bounds(2)])
    ylim ([bounds(3),bounds(4)])
    
    set(gca,'FontSize',15)
    %     set(gcf, 'InvertHardcopy', 'off')
    set(fh,'Units','Inches');
    pos = get(fh,'Position');
    set(fh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fh,'car_sim.pdf','-dpdf','-r0')
end

