classdef two_lane_road < world
    properties
       lane_center_lines
       road_upper_boundary
       road_lower_boundary
       max_obs_speed
       agent_time
       lane_width
       lane_boxes
       road_length
       merge_zone
    end

    methods
function W = two_lane_road(varargin)

    lane_width = 4;
    road_length = 250;
    merge_zone = [150,200];
    verbose = 0;
    bound_space = 2*lane_width;
    
    
    for idx = 1:2:length(varargin)
        switch varargin{idx}
            case 'verbose'
                verbose = varargin{idx+1} ;
            case 'lane_width'
                lane_width = varargin{idx+1} ;
            case 'road_length'
                road_length = varargin{idx+1} ;
            case 'bound_space'
                bound_space = varargin{idx+1} ;
            case 'merge_zone'
                merge_zone = varargin{idx+1};
        end
    end


    bounds = (2*lane_width+bound_space)*[-1 1 -1 1]+[0 road_length 0 0];
    
    W@world('bounds',bounds,'verbose',verbose)
    
    W.lane_center_lines={[0,merge_zone(1);-lane_width/2,-lane_width/2],...
        [0,road_length;lane_width/2,lane_width/2]};
    
    W.lane_width = lane_width;
    
    W.road_length = road_length;
    W.merge_zone = merge_zone;
    merge_length = merge_zone(2)-merge_zone(1);
    
    if merge_zone(2) < road_length
        ypts = linspace(-lane_width,0,ceil(merge_length));
        acoeff = (-lane_width^2+(merge_length))/lane_width;
        xpts =  merge_zone(2)- ypts.*(ypts-acoeff);
     
        W.road_upper_boundary = [0 road_length;lane_width,lane_width];
         W.road_lower_boundary = [0 merge_zone(1) xpts road_length; -lane_width,-lane_width,ypts,0];
    else
        W.road_upper_boundary = [0 road_length;lane_width,lane_width];
        W.road_lower_boundary = [0 road_length;-lane_width,-lane_width];
   
    end
    
    W.start = [0;-lane_width/2];
    W.goal = [road_length;lane_width/2];
    
  
end
    %%
    % obs paths are an Nobs x 1 cell array containing pose positions for
    % the obstacles to move along
function setup(W,obs_paths,varargin)

        W.N_obstacles = length(obs_paths);

        obs_speeds = 20*rand(max([1,W.N_obstacles]))+10;
        refresh = true(W.N_obstacles,1);
        obs_footprint = [4.8,2];
        obs_start_distance = zeros(1,W.N_obstacles);
        obs_start_times = zeros(length(obs_paths),1);

        for idx = 1:2:length(varargin)
            switch varargin{idx}
                case 'obs_speeds'
                    obs_speeds = varargin{idx+1} ;
                case 'refresh'
                    refresh = varargin{idx+1} ;
                case 'obs_footprint'
                    obs_footprint = varargin{idx+1};
                case 'obs_start_distance'
                    obs_start_distance = varargin{idx+1};
            end
        end

        W.max_obs_speed =max(obs_speeds);
        if (W.N_obstacles) == 0
            W.obstacles = [];
        else
            for i=1:W.N_obstacles
                O(i) = railObstacleSimple(obs_footprint(min([i,size(obs_footprint,1)]),:),...
                    obs_paths{i},obs_speeds(min([i,length(obs_speeds)])),...
                    refresh(min([i,length(refresh)])),'start_time',obs_start_times(i));
                
                if obs_start_distance(i)~=0
                    
                    obs_path_length = [0,cumsum(sqrt(diff(obs_paths{i}(1,:)).^2+diff(obs_paths{i}(2,:)).^2))];
                    
                    obs_start_distance(i) = min([obs_start_distance(i),max(obs_path_length)]);
                    obs_start_distance(i) = max([obs_start_distance(i),0]);
                    
                    obs_start_pose = interp_with_angles(obs_path_length',obs_paths{i}',obs_start_distance(i),3)';
                    
                else
                    obs_start_pose = obs_paths{i}(:,1);
                end
                
                O(i).reset(obs_start_pose)
            end

            W.obstacles = O;
        end
      
        W.agent_time =0;    
end

 

%%
function plot_at_time(W,t)
    hold on

    LineWidth=2.0;

    
    %plot boundary
    plot(W.road_lower_boundary(1,:),W.road_lower_boundary(2,:),'k','LineWidth',LineWidth)
    plot(W.road_upper_boundary(1,:),W.road_upper_boundary(2,:),'k','LineWidth',LineWidth)
    
    plot([0 ,W.road_length],[0 0],'--','Color',[0.9290, 0.6940, 0.1250],'LineWidth',LineWidth)
   
    
    for i=1:length(W.obstacles)
        W.obstacles(i).plot_at_time(t);
    end
    
    
    xlim(W.bounds(1:2))
    ylim(W.bounds(3:4))
end
%% reset
function reset(W)
% reset the agent time to 0, and each obstacle's time to 0
reset@world2D(W) ;
W.agent_time = 0 ;
for idx = 1:W.N_obstacles
    W.obstacles(idx).reset() ;
end
end
%% move obstacles
function moveObstacles(W,A)
        % move each obstacle up until the agent's current time
        ta = A.time(end) ;

        for idx = 1:W.N_obstacles
            to = W.obstacles(idx).time(end) ;
            tm = ta - to ;
            if tm > 0
                W.obstacles(idx).move(tm) ;
            end
        end

        % update the world time
        W.agent_time = ta ;
    end

%% crash check
function [out,tcrash] = collision_check(W,A,~)
        W.moveObstacles(A)

        % for each obstacle, find where it is close to the agent's
        % trajectory in both space and time; if there is a crash at these
        % points, check if the agent was stopped
        idx = 1 ;

        % get agent info
        za = A.state ;
        ta = A.time ;
        fa = A.footprint_vertices ;
        dagt = sqrt(sum(A.footprint.^2)) ;

        % create time vector for crash check and reinterpolate agent trajectory
        Nt = ceil((ta(end) - ta(1))/0.005) ;
        tdes = linspace(ta(1),ta(end),Nt) ;
        za = matchTrajectories(tdes,ta,za) ;
       
        
        out = false ;
         tcrash = NaN;
        %check boundary
        dlboundary = dist_points_to_polyline(za(1:2,:),W.road_lower_boundary)<= sqrt(sum(A.footprint.^2));
        zc = za(:,dlboundary);
        tc = tdes(dlboundary);
        idx = 1;
        while out == false && idx <= size(zc,2)
            zfp = zc(A.position_indices,idx)+rotmat(zc(A.heading_index,idx))*A.footprint_vertices;
            xc = polyxpoly(zfp(1,:)',zfp(2,:)',W.road_lower_boundary(1,:)',W.road_lower_boundary(2,:)');
            if ~isempty(xc)
                out = true;
                tcrash = tc(idx);
            end
            idx = idx+1;
        end
        
        duboundary = dist_points_to_polyline(za(1:2,:),W.road_upper_boundary)<= sqrt(sum(A.footprint.^2));
        zc = za(:,duboundary);
        tc = tdes(duboundary);
        idx = 1;
        while out == false && idx <= size(zc,2)
            zfp = zc(A.position_indices,idx)+rotmat(zc(A.heading_index,idx))*A.footprint_vertices;
            xc = polyxpoly(zfp(1,:)',zfp(2,:)',W.road_upper_boundary(1,:)',W.road_upper_boundary(2,:)');
            if ~isempty(xc)
                out = true;
                tcrash = tc(idx);
            end
            idx = idx+1;
        end
        
        
       
        while idx <= W.N_obstacles && out == false
            % get obstacle info
            o = W.obstacles(idx) ;
            zo = o.state;
            to = o.time ;
            fo = o.footprint_vertices ;

            % match the obstacle trajectory to the agent trajectory
            zo = matchTrajectories(tdes,to,zo) ;

            % find where the trajectories are close and the agent is not
            % stopped (meaning traveling < 0.01 cm/s)
            dobs = sqrt(sum(o.footprint.^2)) ;
            
            point_vel=sqrt(diff(zo(o.pose_indices(1),:)).^2+diff(zo(o.pose_indices(2),:)).^2)./diff(tdes);
            
            %logical variabe true if the obstacle is not refreshing at the
            %boundaries of the world
            L_obs_not_refreshing = point_vel<2*abs(o.velocity);
            L_obs_not_refreshing = [L_obs_not_refreshing,L_obs_not_refreshing(end)];
            
            %Local variable true if obstacle is close enough for collision
            %to be possible
            L_close = (sqrt(sum((zo(o.pose_indices([1,2]),:)-za(A.position_indices,:)).^2,1)) < (dobs + dagt));

            dlog = L_close & L_obs_not_refreshing ;
        
            % if any points are potential collisions, place the footprints
            % at those points and check for collisions
            if any(dlog)
                % get the two vectors at the points where they are close
                oc = zo(:,dlog) ;
                ac = za(:,dlog) ;
                tc = tdes(dlog);

                for cidx = 1:size(oc,2)
                    % for each point that may be a crash, place the agent and the
                    % obstacle at that point
                    zoc = oc(o.pose_indices(1:2),cidx) ; % obstacle at potential crash
                    zac = ac(A.position_indices,cidx) ; % agent at potential crash

                    % rotate obstale contour at potential crash point
                    if ~isempty(o.pose_indices(3))
                        roc = oc(o.pose_indices(3),cidx);
                    else
                        roc = 0;
                    end

                    Roc = rotmat(roc);

                    % rotate agent's contour at potential crash point, then
                    % shift contour to point
                    rac = ac(A.heading_index,cidx) ;

                    Rc = rotmat(rac);

                    fac = Rc*fa + repmat(zac,1,size(fa,2)) ;

                    % shift obstacle contour to point
                    foc = Roc*fo + repmat(zoc,1,size(fo,2)) ;

                    % if the two contours intersect, there was a crash!
                    [xc] = polyxpoly(fac(1,:)',fac(2,:)',foc(1,:)',foc(2,:)') ;
                    
                    if ~isempty(xc)
                        
                        %check if agent is going greater than 0.01 m/s
                        L_moving = ac(A.speed_index,cidx)>0.01;
                        
                        if L_moving
                            out = true ;
                            tcrash = tc(cidx);
                        end
                    end
                end
                
            end

            % increment index
            idx = idx + 1 ;
        end
end

%% get nearby obstacles
function O = getNearbyObstacles(W,A,P,~)
        % Method: getNearbyObstacles
        %
        % Return a structure for each obstacle containing a time vector, a
        % position vector, and a footprint

        % prep output
        O = struct([]) ;

        % predict each obstacle's motion
        agent_location = A.state(A.position_indices,end) ;

        tp = P.prediction_horizon ;
        if isempty(tp)
            tp=0;
        end
        Ds = A.sensor_radius ;

        for idx = 1:W.N_obstacles
            % predict future motion of each obstacle
            
            %buffer=P.buffer;
            
            %prediction_time_discretization = P.prediction_time_discretization;
            
            %if isprop(P,'point_spacing')
            %    point_spacing = P.point_spacing ;
            %else
            %    point_spacing = [] ;
            %end
            
            %o = W.obstacles(idx).predict(tp,[],[],[],[],buffer,prediction_time_discretization,point_spacing) ;
            
            %             try
            % if the planner has extra arguments for the predictor send
            % them
            pred_varargin = P.prediction_varargin;
            if isempty(pred_varargin)
                o = W.obstacles(idx).predict(tp,[],[],[],[]) ;
            else
                o = W.obstacles(idx).predict(tp,[],[],[],[],pred_varargin{:}) ;
            end
            
            %catch
            %                o = W.obstacles(idx).predict(tp,[],[],[],[]) ;
            %             end
            
            
            % see if any point of predicted motion is within the sensor
            % horizon of the agent; if so, add the prediction to the output
            % structure (notice that this is actually quite generous to the
            % planner, since the predicted motion could be at the edge of
            % the sensor horizon, but the entire prediction is passed to
            % the planner)
            Zpred = o.state(W.obstacles(idx).pose_indices(1:2),:) ;
            d = distPointToPoints(agent_location,Zpred) ;
            if any(d <= Ds)
                O = [O;o] ;
            end
        end
        
  
    end

%% get workd info
function W_info = get_world_info(W,agent_info,P)
    W_info.obstacle_predictions = W.getNearbyObstacles(agent_info,P);
    W_info.bounds = W.bounds;
    W_info.lane_boxes = W.lane_boxes;
    W_info.max_obs_speed = W.max_obs_speed;
    W_info.lane_width = W.lane_width;
    W_info.road_upper_boundary = W.road_upper_boundary;
    W_info.road_lower_boundary = W.road_lower_boundary;
end
    end

end
