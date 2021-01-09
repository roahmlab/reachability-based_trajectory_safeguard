function [dmin,dout,pout,pidx] = distPointToPolyline(p,P)
% Function: [dmin,dout,pout,pidx] = distPointToPolyline(p,P)
%
% Given a point p (2-by-1) and a polyline (2-by-m), return a 1-by-m vector
% of the closest distances from the point to the polyline and (optionally)
% the minimum distance, the closest point on the polyline to p, and the
% index of the line segment on which the closest point lies (if the index
% is m, then the closest point is the last point of the polyline P)
%
% The polyline is P consists of m-1 line segments, from the first to m-1
% points "Pa" to the second to mth points "Pb" so the output distances are
% to these m-1 line segments and to the last point of Pb.
%
% Modified from: http://www.alecjacobson.com/weblog/?p=1486

    Pa = P(:,1:end-1) ;
    Pb = P(:,2:end) ;

    dP = Pb - Pa ;
    P2 = sum(dP.*dP,1) ;
    
    P2p = repmat(p,1,size(P,2)-1) - Pa ;

    t = sum(P2p.*dP,1)./P2 ;

    tlog = t > 0 & t < 1 ;

    if any(tlog)
        Pa(:,tlog) = Pa(:,tlog) + repmat(t(tlog),2,1).*dP(:,tlog) ;
        Pall = [Pa,P(:,end)] ;
    else
        Pall = P ;
    end
    
    dout = distPointToPoints(p,Pall) ;
    [dmin,pidx] = min(dout) ;
    
    if nargout > 2
        pout = Pall(:,pidx) ;
    end
end