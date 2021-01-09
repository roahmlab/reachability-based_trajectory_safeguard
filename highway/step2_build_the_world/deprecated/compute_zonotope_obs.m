function [O_local] = compute_turtlebot_discretized_obs(O_world,turtlebot_pose,b,r,FRS,O_bounds)
% O_FRS = compute_turtlebot_discretized_obs(O_world,turtlebot_pose,b,r,FRS)
% O_FRS = compute_turtlebot_discretized_obs(O_world,turtlebot_pose,b,r,FRS,O_bounds)
%
% Take obstacle points defined in a global coordinate frame and transform
% them into the scaled, shifted FRS frame for the TurtleBot.
%
% Author: Shreyas Kousik
% Created: 30 May 2019

    % for now, just don't use the arc point spacing by setting the miter limit
    % to 2 in buffer_polygon_obstacles
    O_buf = buffer_polygon_obstacles(O_world,b,2) ;
    
    % handle the world bounds obstacle if it's passed in
%     if nargin > 5
%         O_bounds_buf = buffer_polygon_obstacles(O_bounds,b,2) ;
%         O_buf = [O_buf, nan(2,1), O_bounds_buf] ;
%     end
%     
%     % create discretized obstacle points
%     O_pts = interpolate_polyline_with_spacing(O_buf,r) ;
% plot(O_buf(1,:),O_buf(2,:))
%     P.current_obstacles
    % put the obstacles into the FRS frame
%     x0 = FRS.initial_x ;
%     y0 = FRS.initial_y ;
%     D = FRS.distance_scale ;
%     O_FRS = world_to_FRS(O_pts,turtlebot_pose,x0,y0,D) ;
    O_local = world_to_local(turtlebot_pose,O_buf);
    % filter out points that are too far away to be reached
%     O_FRS = crop_points_outside_region(0,0,O_FRS,1) ;
    O_local = O_local';
    % return the buffered obstacle as well
end