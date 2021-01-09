%% get Euclidean distances between points along a path
function [cumd,d] = cumulativeDist(p)
    % For n points in p (n x m), get the distances from p_i to p_{i+1} and
    % return d as an n x 1 vector
    d = diff(p,1) ;
    d = [0;sqrt(sum(d.*d,2))] ;
    cumd = cumsum(d) ;
end