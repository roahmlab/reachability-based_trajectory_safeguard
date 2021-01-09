function [Xo,Yo] = shiftScaleAndRotateContourAxes(X,Y,z,h,D,z0)
    if nargin < 6
        z0 = [0;0] ;
    end
    
    x0 = z0(1) ; y0 = z0(2) ;
    
    Xo = D.*((X - x0/D).*cos(h) - (Y - y0/D).*sin(h)) + z(1) ;
    Yo = D.*((X - x0/D).*sin(h) + (Y - y0/D).*cos(h)) + z(2) ;
end