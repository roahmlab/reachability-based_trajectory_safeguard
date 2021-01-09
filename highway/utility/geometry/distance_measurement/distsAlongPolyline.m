function dP = distsAlongPolyline(P,first_value)
% Given a polyline P (n-by-N) in n dimensions of N points, return a 1-by-N
% vector dP in which the first entry is 0 (by default, or can be given by
% the user), and each ith entry is the distance from the (i-1)th point to
% the ith point.
    if nargin < 2
        first_value = 0 ;
    end
    
    if ~isempty(first_value)
        dP = [first_value,sqrt(sum(P.^2,1))] ;
    else
        dP = sqrt(sum(P.^2,1)) ;
    end
end