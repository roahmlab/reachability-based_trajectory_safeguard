function [O,Plog] = filterPoints(x,y,O,L,method)
    % Given a n-by-n set of points, a position (x,y), and a distance L, remove
    % all points that are not within the [-L,L]^2 box, or not in the L-radius
    % super-ellipse. The 5th input parameter can be either 'rect' or an even
    % number specifying the degree of the super-ellipse
    
    if isempty(O)
        O = [] ;
        Plog = [] ;
    else
        P = O(1:2,:) ;

        if nargin < 5
            method = 'rect' ;
        end

        if strcmp(method,'rect')
            Plog = any(abs(P - repmat([x;y],1,size(P,2))) > L, 1) ;
        else % use a super-ellipse
            d = sqrt(sum((abs(P-repmat([x;y],1,size(P,2)))).^method,1)) ;
            Plog = d > L ;
        end

        O = O(:, ~Plog) ;
    end
    
    if nargout < 2
        Plog = [] ;
    end
end
