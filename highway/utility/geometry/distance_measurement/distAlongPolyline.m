function [d,pct] = distAlongPolyline(p,P)
% Given a point p (2-by-1) on a polyline P (2-by-n), return the distance of
% the point along the polyline, and (optionally) the percentage of the
% total distance
    [~,~,p,pidx] = distPointToPolyline(p,P) ;
    pct=[];
    
    % get distances along the polyline
    Pa = P(:,1:end-1) ; Pb = P(:,2:end) ;
    dP = Pb - Pa ;
    dPdists = sqrt(sum(dP.^2,1)) ;
    Pdists = [0, cumsum(dPdists)] ;
    Ppcnts = Pdists./Pdists(end) ;

    % get distance up to closest index and add distance from the point at
    % the closest index to the point itself
    d0 = Pdists(pidx) ;
    d1 = distPointToPoints(p,P(:,pidx)) ;
    d = d0 + d1 ;
    
    if nargout > 1
        p0 = Ppcnts(pidx) ;    
        p1 = d1./Pdists(end) ;
        pct = p0 + p1 ;
    end
end