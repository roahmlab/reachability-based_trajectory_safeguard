function dout = distPointsToPolyline(points,polyline)
% Function: dout = distPointsToPolyline(points,polyline)
%
% Given points p (2-by-n) and a polyline (2-by-m), return a 1-by-m vector
% of the closest distances from the points to the polyline.
%
% The polyline is P consists of m-1 line segments, from the first to m-1
% points "Pa" to the second to mth points "Pb" so the output distances are
% to these m-1 line segments and to the last point of Pb.
%
% Modified from: http://www.alecjacobson.com/weblog/?p=1486

    dout = arrayfun(@(k) distPointToPolyline(points(:,k),polyline), 1:size(points,2)) ;
end