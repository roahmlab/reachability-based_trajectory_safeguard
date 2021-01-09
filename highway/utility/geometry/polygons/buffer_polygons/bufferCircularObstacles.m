function O = bufferCircularObstacles(c,r,b,l)
% Given a list of centers (2-by-N) and radii (1-by-N), and a buffer
% distance b (scalar), create polylines that buffer each circle by the
% distance b, with points on the polyline no further than the distance l
% apart

O = [] ;

for idx = 1:size(c,2)
    ci = c(:,idx) ;
    ri = r(:,idx) + b;
    
    % determine the number of points required given the current circle's
    % size
    N = ceil(2*pi*ri/l) ;
    t = linspace(0,2*pi,N) ;
    
    % create the required polyline
    o = [r.*cos(t) + ci(1) ;
         r.*sin(t) + ci(2)] ;
     
    O = [O, nan(2,1), o] ;    
end

% move initial column of nans to end
O = [O(:,2:end), nan(2,1)] ;
end