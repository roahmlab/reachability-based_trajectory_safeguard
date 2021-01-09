function h = headingPoseToPoints(p,P)
% h = headingPointToPoints(p,P)
%
% Given a pose (x,y,h) in SE^2 and points P (2-by-N), return a 1-by-N
% vector of the heading in (-pi,pi] from p to each of the points in P
    h = headingPointToPoints(p(1:2),P) - p(3) ;
end