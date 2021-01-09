function P = loc2glob(Q,z,h,D,x0,y0)
    % Method: q = loc2glob(Q,z,h,D,x0,y0)
    %
    % Rotate and shift the local points Q (2-by-n) to the global frame
    % of the segway, given the current position z and heading h

    % make rotation matrix
    R = [cos(h), sin(h) ; -sin(h), cos(h)] ;

    % create output
    N = size(Q,2) ;
    P = D.*R'*(Q - repmat((1/D).*[x0;y0],1,N)) + repmat(z,1,N) ;
end