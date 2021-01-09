function Safe=hit_obs(state,Obs,Zono , plot_flag,color)

[xi yi]= polyxpoly(Obs(:,1),Obs(:,2),Zono(:,1),Zono(:,2));
if plot_flag
 figure(1)
 Obs = local_to_world(state,Obs')';
 Zono = local_to_world(state,Zono')';
%  hobs = mapshow(Obs(:,1),Obs(:,2),'DisplayType','polygon','FaceAlpha',1,'LineStyle',':');
 if ~isempty(xi)
     color = [1 0 0];
 end
 if color(1)<0.3
     hzono= fill(Zono(:,1),Zono(:,2),color);
      hzono.FaceAlpha = 0.4;
      hzono.EdgeAlpha = 0;
%       uistack(hzono,'bottom')
%  hzono= mapshow(Zono(:,1),Zono(:,2),'DisplayType','line', 'LineStyle','-','Color',color);
 else
      hzono= mapshow(Zono(:,1),Zono(:,2),'DisplayType','line','Color',color);
      
%       hzono.
 end
%  mapshow(xi,yi,'DisplayType','point','Marker','o')
% drawnow
%     pause(0.01);
%     delete(hobs);
 %    delete(hzono);

end

if ~ isempty(xi) 
    Safe = false;
else
    Safe = true;
end 
end