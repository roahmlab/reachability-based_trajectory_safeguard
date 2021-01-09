function P = makeRandomPolygon(Nvertices,offset,scale)
% P = makeRandomPolygon(Nvertices,offset,scale)
%
% Creates a random 2-D polygon with a given number of vertices, offset from
% the point (0,0), and scaled up or down as desired. The polygon is a
% 2-by-N, counterclockwise polyline.
%
% If the inputs Nvertices, offset, and scale are vectors, then the output
% is a 2-by-N polyline with as many polygons are the length of Nvertices,
% with each polygon separated by a column of NaNs.

    if nargin < 1
            Nvertices = round(randRange(3,9)) ;
    end
    
    Npolygons = size(Nvertices,2) ;
    
    if nargin < 2
        offset = zeros(2,Npolygons) ;
    end
    
    if nargin < 3
        scale = ones(1,Npolygons) ;
    end

    if length(Nvertices) == 1
        vertices = rand(2,Nvertices) - 0.5 ;
        P = scale*pointsToCCW(vertices, true) + offset ;
    else
        P = [] ;
        for idx = 1:length(Nvertices)
            P = [P, nan(2,1), makeRandomPolygon(Nvertices(idx),offset(:,idx),scale(idx))] ;
        end
        P = P(:,2:end) ;
    end
end