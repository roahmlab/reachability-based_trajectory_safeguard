function in = inpoly_multipoly(P,O)
% in = inpoly_multipoly(P,O)
%
% Run inpoly on multiple polygons. The input P is a 2-by-n vector of (x,y)
% points to test for being in the polygons in O, which is a 2-by-n vector
% of CCW polygons separated by NaNs. This function assumes there are no
% polygons with holes in them (which would be CW).

    % get x and y points and polygons
    Px = P(1,:)' ;
    Py = P(2,:)' ;
    Ox = O(1,:)' ;
    Oy = O(2,:)' ;
    
    % pad O with nans if needed
    if ~isnan(O(1,end))
        O = [O, nan(2,1)] ;
    end
    
    % get indices of nans, i.e. of separate polygons
    nanlog = isnan(Ox) ;
    nanidx = find(nanlog) ;
    idxs = [nanidx-1,nanidx+1]' ;
    idxs = [1; idxs(:)] ;
    
    % run insidepoly on each polygon
    in = false(size(Px)) ;
    for idx = 1:2:(length(idxs)-1)
        ii = idxs(idx:idx+1) ;
        otemp = O(:,ii(1):ii(2)) ;
        intemp = insidepoly(Px,Py,otemp(1,:)',otemp(2,:)') ;
        in = in | intemp ;
    end
end