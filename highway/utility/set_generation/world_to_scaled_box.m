function var_box = world_to_scaled_box(var_world,world_range,box_range)
N = length(var_world);
if nargin<3
    box_range = [-ones(N,1) ones(N,1)];
end
if numel(box_range) ==1
    box_range = abs(box_range)*[-ones(N,1) ones(N,1)];
end

var_box = (1./((world_range(:,2)-world_range(:,1)))) .* (var_world-world_range(:,1)).*(box_range(:,2)-box_range(:,1)) + box_range(:,1);


end

