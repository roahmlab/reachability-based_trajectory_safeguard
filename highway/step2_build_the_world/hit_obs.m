function Safe=hit_obs(state,Obs,Zono , plot_flag,color)
% This is a helper function to check for obstacle intersection naively!
% it is also helpeful to plot the zonotopes to visualize bugs.
[xi yi]= polyxpoly(Obs(:,1),Obs(:,2),Zono(:,1),Zono(:,2));
if plot_flag
 figure(1)
 Obs = local_to_world(state,Obs')';
 Zono = local_to_world(state,Zono')';

 if ~isempty(xi)
     color = [1 0 0];
 end
 if color(1)<0.3
     hzono= fill(Zono(:,1),Zono(:,2),color);
      hzono.FaceAlpha = 0.4;
      hzono.EdgeAlpha = 0;

 else
      hzono= mapshow(Zono(:,1),Zono(:,2),'DisplayType','line','Color',color);
      
%       hzono.
 end

end

if ~ isempty(xi) 
    Safe = false;
else
    Safe = true;
end 
end