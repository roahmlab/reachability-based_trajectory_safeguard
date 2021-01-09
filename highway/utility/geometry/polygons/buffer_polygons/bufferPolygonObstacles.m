function Pout = bufferPolygonObstacles(P,b,miterlim)
% Given obstacles P as a 2-by-N array of (x,y) points, buffer them by a
% distance b and return a 2-by-(however many are needed) polygon array

  P = polyshape(P(1,:),P(2,:)) ;
    if nargin <3
        P = polybuffer(P,b) ;   
    elseif strcmp(miterlim,'square')
        P = polybuffer(P,b,'JointType','square');
    else
        P = polybuffer(P,b,'JointType','miter','MiterLimit',miterlim);
    end
  Pout = P.Vertices' ;
end