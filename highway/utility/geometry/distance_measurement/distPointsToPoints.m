function [D,d] = distPointsToPoints(P1,P2)
    % given two sets of points P1 (2 x m) and P2 (2 x n), calculate the
    % distances from every point of P1 to every point of P2 and return it
    % as a matrix D (m x n) and as a vector (1 x mn)
    [D,d] = distManyToMany(P1,P2) ;
end