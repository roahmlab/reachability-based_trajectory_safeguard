function points_out = rotatePointsToZeroHeading(pose,points_in)
% points_out = rotatePointsToZeroHeading(pose,points_in)
%
% Given a pose in SE^2 as (x,y,h), transform the points passed in so that
% they are "centered" at the provided pose and rotated so that the angle of
% the provided pose is set to 0

    R = rotmat(-pose(3)) ;
    Z = pose(1:2) ;
    
    points_out = R*(points_in - Z) + Z ;
end