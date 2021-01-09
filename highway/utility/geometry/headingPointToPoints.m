function h = headingPointToPoints(p,P)
% h = headingPointToPoints(p,P)
%
% Given a point p (2-by-1) and points P (2-by-N), return a 1-by-N vector of
% the heading in (-pi,pi] from p to each of the points in P
    h = atan2(P(2,:) - p(2), P(1,:) - p(1)) ;
end