function [discrete_points,buffered_polygon_vertices] = discretize_polygon_obstacles(P_vertices,buffer,boundary_point_spacing,interior_point_spacing)


if nargin < 4
    interior_point_spacing = Inf;
    if nargin < 3
        boundary_point_spacing = Inf;
        if nargin < 2
            buffer = 0;
        end
    end
end

if isempty(buffer)
    buffer = 0;
end
if isempty(boundary_point_spacing)
    boundary_point_spacing = Inf;
end
if isempty(interior_point_spacing)
    interior_point_spacing = Inf;
end


buffered_polygon_vertices = bufferPolygonObstacles(P_vertices,buffer,'square');

nan_idxs = find(isnan(buffered_polygon_vertices(1,:)));

nan_idxs = [0,nan_idxs,size(buffered_polygon_vertices,2)+1];

discrete_points = cell(1,length(nan_idxs));

for i = 1:length(nan_idxs)-1
    
    boundary_vertices = buffered_polygon_vertices(:,(nan_idxs(i)+1) : (nan_idxs(i+1)-1));
    
    %increase density according to point spacing
    if  boundary_point_spacing < Inf
        %close polygon for point spacing algorithm to work
        if norm(boundary_vertices(:,end)-boundary_vertices(:,1))>0
            boundary_vertices=[boundary_vertices,boundary_vertices(:,1)];
        end
        
        poly_points = increasePolylineDensity(boundary_vertices,boundary_point_spacing);
        
    else
        poly_points = boundary_vertices;
    end
    
    %add points for interior
    if interior_point_spacing < Inf
        
        grid_points_to_add = hex_grid(poly_points(1,:),poly_points(2,:),interior_point_spacing);
        
        if ~isempty(grid_points_to_add)
            
            L_pts_in_polygon = inpolygon(grid_points_to_add(1,:)',grid_points_to_add(2,:)',boundary_vertices(1,:)',boundary_vertices(2,:)')';
            
            pts_to_add_to_slice = grid_points_to_add(:,L_pts_in_polygon );
            
            poly_points = [poly_points,pts_to_add_to_slice];
        end
        
    end
    
    %unclose and elimiate anypoints that are too close together
    poly_points = uniquetol(poly_points','ByRows',true)';
    
    discrete_points{i} = poly_points;
end

discrete_points = cat(2,discrete_points{:});


end