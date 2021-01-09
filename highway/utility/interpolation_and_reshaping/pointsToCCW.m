function Pout = pointsToCCW(Pin,close_loop)
% Given a 2-by-n set of points where the first row is x coordinates and the
% second row is y coordinates, rearrange the points in counter-clockwise
% order. If the "close loop" flag is set to true
    if nargin < 2
        close_loop = false ;
    end
    
    % get x/y coords
    Px = Pin(1,:) ; Py = Pin(2,:) ;
    
    % get centroid
    cx = mean(Px) ; cy = mean(Py) ;
    
    % get angle and distance from centroid to each point
    ang = atan2(Py-cy,Px-cx) ;
    dst = distPointToPoints([cx;cy],Pin) ;
    
    % sort by angle in [-pi,pi] then by distance
    [~,idx] = sortrows([ang(:) dst(:)]) ;
    
    % rearrange points and create output
    Pout = [Px(1,idx) ; Py(1,idx)] ;
    if close_loop
        Pout = [Pout, Pout(:,1)] ;
    end
end