function q = glob2loc(P,z,h,D,x0,y0)
    % Method: q = glob2loc(P,z,h)
    %
    % Rotate and shift the global points P (2-by-n) to the local frame
    % of the segway, given the current position z and heading h
    
    % make rotation matrix
    R = [cos(h), sin(h) ; -sin(h), cos(h)] ;

    % create output
    q = R*((1/D)*(P - repmat(z,1,size(P,2)))) + repmat((1/D)*[x0;y0],1,size(P,2)) ;
end