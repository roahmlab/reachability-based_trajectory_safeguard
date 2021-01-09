function out = isVertex(v,V,posflag)
% Given a point v (m-by-1) and a collection of points V (m-by-n), return
% true if v is a column of V or false otherwise. If the third argument is
% set to true, then the output is a 1-by-n logical that is true for each
% column of V that is given by v.
    if nargin < 3
        posflag = false ;
    end
    
    N = size(V,2) ;
    vmat = repmat(v,1,N) ;
    
    if posflag
        out = all(vmat == V,1) ;
    else
        out = any(all(vmat == V,1)) ;
    end
end