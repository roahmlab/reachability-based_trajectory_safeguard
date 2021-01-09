function Pout = getNearbyPolygons(p,P,sensor_region)
% Pout = getNearbyPolygons(p,P,distance)
% Pout = getNearbyPolygons(p,P,zone)
%
% Given a point p (x,y), polygons P (2-by-N points separated by NaNs), and
% a distance d, return all polygons with at least one point within a
% distance d of p.
%
% If the last argument is a polygonal (CCW) region, use that instead of the
% distance d. Note that this does not use p (the first input).

    if isscalar(sensor_region)
        dlog = (distPointToPoints(p(1:2),P) <= sensor_region) ;
    else
        [in,on] = inpolygon(P(1,:)',P(2,:)',...
                            sensor_region(1,:)',sensor_region(2,:)') ;
        dlog = in' | on' ;
    end

    dlog = dlog | isnan(P(1,:)) ;
    
    Pout = P(:,dlog) ;
end