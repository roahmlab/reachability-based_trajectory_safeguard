% create agent and world
A = quadrotor_agent() ;
W = zonotope_box_world_safe_rl('start',[50;2;5],'buffer_start',0) ;

A.reset(W.start) ;

ob = W.get_ob(A.get_agent_info()) ;

% get obstacle distance observation data
D = ob(10:end) ;
D_max = max(D) ;
D_min = min(D) ;

% get P_ob
P_ob = W.obstacle_observation_points + repmat(W.start,1,length(D)) ;

%% plotting
figure(1) ; clf ; hold on ; axis equal ; view(3) ;

plot(W) ;
plot_path(W.start,'rp')

for idx = 1:length(D)
    r = (1-(D(idx) - D_min)./(D_max - D_min)).*[1 0 0] ;
    g = (D(idx) - D_min)./(D_max - D_min).*[0 1 0] ;
    plot_path(P_ob(:,idx),'.','color',r + g) ;
end