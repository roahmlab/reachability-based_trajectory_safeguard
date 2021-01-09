function zono = local_to_zono(vertices)
x_g = 0.5*(max(vertices(:,1)) - min(vertices(:,1)));
x_center = 0.5*(max(vertices(:,1)) + min(vertices(:,1)));
y_g = 0.5*(max(vertices(:,2)) - min(vertices(:,2)));
y_center = 0.5*(max(vertices(:,2)) + min(vertices(:,2)));

obs_center = [x_center; y_center];
obs_gen = [[x_g; 0], [0; y_g]];
zono = zonotope([obs_center, obs_gen]);
end