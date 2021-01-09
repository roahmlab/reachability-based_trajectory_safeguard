function plot_agent_pose(A,Z,idxs,cs,body_vargin,arrowcs,arrow_vargin)
if nargin<6
    if nargin<5
        body_vargin={'EdgeColor','k'};
        if nargin<4
            cs = [0.8,0.8,1.];
        end
    end
    arrow_vargin = body_vargin;
    arrowcs = cs;
end
for i=idxs
   if isempty(A.heading_state_index)
       h=0;
   else
    h=Z(A.heading_state_index,i);
   end
    R=rotmat(h);
    fc=R*A.footprint_contour+[Z(A.xy_state_indices(1),i);Z(A.xy_state_indices(2),i)];
    hold on
    fill(fc(1,:),fc(2,:),cs,body_vargin{:})
    if ~isempty(A.heading_state_index)
        thetas = linspace(0,2*pi,4) ;
        
        if length(A.footprint) == 2
            D = max(A.footprint) / 2 ;
            trans = mean([max(A.footprint_contour(1,:)),min(A.footprint_contour(1,:))]);
            xdm=1/3;
        else
            D = A.footprint ;
            trans=0;
            xdm=1/2;
        end
        
        arrow = R*[((D*xdm).*cos(thetas) +trans) ;...
            (D/4).*sin(thetas)]+ [Z(A.xy_state_indices(1),i);Z(A.xy_state_indices(2),i)];
        hold on
        fill(arrow(1,:),arrow(2,:),arrowcs,arrow_vargin{:})
    end
end


end 