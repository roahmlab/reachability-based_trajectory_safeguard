classdef highway_HLP < high_level_planner
methods
    function HLP = highway_HLP(varargin)
        HLP@high_level_planner(varargin{:}) ;
    end

    function waypoint = get_waypoint(HLP,agent_info,world_info,d)
        % d is a 3x1 matrix that holds the front car info.
        [~, laneidx] = min(d);
%         [distf, ~   ] = min(d);
        dt = 2; v_max = 4;
        ypos = agent_info.state(2,end); % get lane with furthrest car
        if laneidx == 1
            desired_y = 10;
        elseif laneidx == 2
            if ypos > 6
                desired_y = 10;
            else
                desired_y = 2;
            end
        elseif laneidx == 3
            desired_y = 2;
        end
        
%         desired_y = 4*(laneidx-1)+2;
        desired_x = dt*v_max+agent_info.state(1,end);%*double(distf>30)+dt*v_max/2*double(distf<=30)+
        waypoint = [desired_x;desired_y];

        % update current waypoints
        HLP.current_waypoint = waypoint ;
        HLP.waypoints = [HLP.waypoints, waypoint] ;
        HLP.current_waypoint_index = HLP.current_waypoint_index + 1 ;
    end
end
end