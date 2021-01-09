function [x,y] = rotxy(xin,yin,h)
    R = rotmat(h) ;
    Z = R*[xin;yin] ;
    x = Z(1,:) ;
    y = Z(2,:) ;
end