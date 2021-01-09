function P_out = world_2_local_2scale(P_robot, P_world, set_dx, set_dy, Dx,Dy)
% P_out = world_2_local_2scale(P_robot, P_world, set_dx, set_dy, Dx,Dy)
%
% Given a position and heading of the robot in the world frame, positions
% of obstacles as xy points, and an offset of the robot in its own local
% frame (set_dx, set_dy), transform the points from the world frame to the
% shifted local frame and scale them down by the distance D
%
% INPUTS
%   P_robot     robot position (x,y,heading)
%   P_world     obstacle points in world (x,y) 2 x Nobs
%   set_dx      x position of robot in it's local frame
%   set_dy      y position of robot in it's local frame
%   Dx          scaling in the x dimension
%   Dy          scaling in the y dimension (equal to Dx by default)
%
% OUTPUTS:
%   P_out       obstacle points in robot's local frame

    if ~exist('Dy','var')
        Dy = Dx ;
    end
    
    % extract position and heading from input
    X = P_robot(1,1) ;
    Y = P_robot(2,1) ;
    h = P_robot(3,1) ;
    
    % get the number of world points
    [N_rows, N_cols] = size(P_world);
    
    % shift all the world points to the position of the robot and scale
    % them down to the provided distance
    Imat = ones(N_rows,N_cols);
    
    XY = [X, 0;
          0, Y];
      
    P_out = (P_world - XY*Imat);
    
    % create the rotation matrix to use on these shifted points
    R = [cos(h), sin(h);
         -sin(h), cos(h)];
     
    % create the final shift to move things in local coordinates after
    % rotation
    dXY = [set_dx, 0;
           0, set_dy];
    
    % create the final shifted and scaled version of the obstacle points
    P_out = dXY*Imat + R*P_out;
    
    P_out = [(1/Dx)*P_out(1,:) ; 
             (1/Dy)*P_out(2,:)] ;
 end
