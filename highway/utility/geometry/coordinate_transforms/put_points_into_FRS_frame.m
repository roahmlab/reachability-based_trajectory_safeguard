function P_out = put_points_into_FRS_frame(P_robot, P_world, FRS)
% P_out = put_points_into_FRS_frame(P_robot, P_world, FRS)
%
% Given a position and heading of the robot in the world frame, points or
% poses as (x,y) or (x,y,h), put those points/poses into the FRS frame
%
% INPUTS
%   P_robot     robot position (x,y,heading)
%   P_world     obstacle points in world (x,y) 2 x Nobs
%   FRS         FRS information (see below)
%
% INPUTS EXTRACTED FROM FRS
%   x0          x position of robot in it's local frame
%   y0          y position of robot in it's local frame
%   Dx          scaling in the x dimension
%   Dy          scaling in the y dimension (equal to Dx by default)
%
% OUTPUTS:
%   P_out       points in FRS frame

    % extract position and heading from input
    x = P_robot(1,1) ;
    y = P_robot(2,1) ;
    h = P_robot(3,1) ;
    
    % get the number of world points
    [N_rows, N_cols] = size(P_world);
    
    % shift all the world points to the position of the robot and scale
    % them down to the provided distance
    I_mat = ones(2,N_cols);
    
    xy_mat = [x, 0;
              0, y];
      
    P_out = (P_world(1:2,:) - xy_mat*I_mat);
    
    % create the rotation matrix to use on these shifted points
    R = rotmat(-h) ;
    
    % create the final shift to move things in local coordinates after
    % rotation
    dXY = [FRS.x0, 0;
           0, FRS.y0];
    
    % create the final shifted and scaled version of the obstacle points
    P_out = dXY*I_mat + R*P_out;
    
    P_out = [(1/FRS.Dx)*P_out(1,:) ; 
             (1/FRS.Dy)*P_out(2,:)] ;
         
    % add heading if P_world is poses instead of points
    if N_rows == 3
        P_out = [P_out ; P_world(3,:) - h] ;
    end
 end
