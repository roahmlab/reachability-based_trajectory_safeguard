function plot_agent_and_obstacles_at_time_range(Agent,Obstacles,t_range,figure_num,P)

if nargin<4
    figure
    hold on
else
    figure(figure_num)
    hold on
 
end

%% plot agents

Z_agent = interp1(Agent.time',Agent.state',t_range(:),'pchip')';


plot_agent_pose(Agent,Z_agent,1:size(Z_agent,2),'b')

for i=1:length(Obstacles)
    
    Z_obs =  interp1(Obstacles(i).time',Obstacles(i).state',t_range(:),'pchip')';
    
    
    
    plot_agent_pose(Obstacles(i),Z_obs,1:size(Z_obs,2),'r')

end

%% if we have obstacle_points plot planner obs pts closest to time
if nargin >=5

for i=1:length(t_range)
%find iteration for current planning time
iter_idx = find((t_range(i)-P.info.agent_time)>=0,1,'last');
FRS_idx = P.info.FRS_index_used_for_plan(iter_idx);
Tscale = P.FRS{FRS_idx}.Tscale;
Toffset = P.FRS{FRS_idx}.Toffset;
current_obstacles = P.info.current_obstacles{iter_idx};

%find  abd plot obstacle slices closest to t_range
idx_before = find(t_range(i)-P.info.agent_time(iter_idx)-current_obstacles(3,:)>=0,1,'last');

t_before = current_obstacles(3,idx_before);

Lbefore=abs(current_obstacles(3,:) - t_before)<1e-9;

plot(current_obstacles(1,Lbefore),current_obstacles(2,Lbefore),'.','Color',[1,0.7,0.9])

idx_after = find(t_range(i)-P.info.agent_time(iter_idx)-current_obstacles(3,:)<0,1,'first');

if~isempty(idx_after)
    
    t_after = current_obstacles(3,idx_after);

    Lafter=abs(current_obstacles(3,:) - t_after)<1e-9;

    plot(current_obstacles(1,Lafter),current_obstacles(2,Lafter),'.','Color',[1,0.7,0.9])
end

% plot contours at those times
for k=1:length(Obstacles)
    
    Z_obs =  interp1(Obstacles(k).time',Obstacles(k).state',[t_before;t_after]+P.info.agent_time(iter_idx),'pchip')';
    
    
    
    plot_agent_pose(Obstacles(i),Z_obs,1:size(Z_obs,2),'r')

end

Z_agent = interp1(Agent.time',Agent.state',[t_before;t_after]+P.info.agent_time(iter_idx),'pchip')';


plot_agent_pose(Agent,Z_agent,1:size(Z_agent,2),'b')


for j =1:length(P.FRS{FRS_idx}.out)
%     
    t_scaled = (t_range(i)-P.info.agent_time(iter_idx)+Toffset(j))/Tscale(j);
    %plot frs_contours
    if t_scaled >=0 && t_scaled <=1
        wk = msubs(P.FRS{FRS_idx}.out(j).w,[P.FRS{FRS_idx}.t{j};P.FRS{FRS_idx}.k{j}],[t_scaled;P.info.k_opt_found(:,iter_idx)]);
        
        plot_w_contour(wk,P.FRS{FRS_idx}.z{j},P.info.agent_pose_used_for_plan(:,iter_idx),...
            P.FRS{FRS_idx}.x0(j),P.FRS{FRS_idx}.y0(j),P.FRS{FRS_idx}.Dx(j),P.FRS{FRS_idx}.Dy(j),[0.240000 0.700000 0.440000],'-')
    end  
    
    t_before_scaled = (t_before+Toffset(j))/Tscale(j);
    %plot frs_contours
    if t_before_scaled>=0 && t_before_scaled <=1
        wk = msubs(P.FRS{FRS_idx}.out(j).w,[P.FRS{FRS_idx}.t{j};P.FRS{FRS_idx}.k{j}],[ t_before_scaled;P.info.k_opt_found(:,iter_idx)]);
        
        plot_w_contour(wk,P.FRS{FRS_idx}.z{j},P.info.agent_pose_used_for_plan(:,iter_idx),...
            P.FRS{FRS_idx}.x0(j),P.FRS{FRS_idx}.y0(j),P.FRS{FRS_idx}.Dx(j),P.FRS{FRS_idx}.Dy(j),[0.240000 0.700000 0.440000],':')
    end 
    
    t_after_scaled = (t_after+Toffset(j))/Tscale(j);
    %plot frs_contours
    if t_after_scaled>=0 && t_after_scaled <=1
        wk = msubs(P.FRS{FRS_idx}.out(j).w,[P.FRS{FRS_idx}.t{j};P.FRS{FRS_idx}.k{j}],[ t_after_scaled;P.info.k_opt_found(:,iter_idx)]);
        
        plot_w_contour(wk,P.FRS{FRS_idx}.z{j},P.info.agent_pose_used_for_plan(:,iter_idx),...
            P.FRS{FRS_idx}.x0(j),P.FRS{FRS_idx}.y0(j),P.FRS{FRS_idx}.Dx(j),P.FRS{FRS_idx}.Dy(j),[0.240000 0.700000 0.440000],':')
    end 
    
    
end
    
end
end

end