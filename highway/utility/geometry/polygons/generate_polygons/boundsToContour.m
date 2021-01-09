function C = boundsToContour(B)
    xlo = B(1) ; xhi = B(2) ; ylo = B(3) ; yhi = B(4) ;
    C = [xlo xhi xhi xlo xlo ; ylo ylo yhi yhi ylo] ;
end
