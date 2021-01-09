function plot_ref3(AH,ref, proposed_ref)
Z_prime = ref;
if ~isempty(proposed_ref)
    Z = proposed_ref ;
    
    % plot current plan
    if check_if_plot_is_available(AH,'proposed_trajectory')
        AH.plot_data.proposed_trajectory.XData = Z(1,:) ;
        AH.plot_data.proposed_trajectory.YData = Z(2,:) ;
        AH.plot_data.proposed_trajectory.ZData = Z(3,:) ;

        AH.plot_data.proposed_traj2.XData = Z(1,:) ;
        AH.plot_data.proposed_traj2.YData = Z(2,:) ;
        AH.plot_data.proposed_traj2.ZData = Z(3,:) ;

        %             AH.plot_data.proposed_trajectory.ZData = Z(3,:) ;
    else
        trajectory_data = plot3(Z(1,:),Z(2,:),Z(3,:),'k-','LineWidth',3);
        traj2           =plot3(Z(1,:),Z(2,:),Z(3,:),'Color','y','LineWidth',3,'LineStyle','--');
        AH.plot_data.proposed_trajectory = trajectory_data ;
        AH.plot_data.proposed_traj2 = traj2;
    end
end
Z = Z_prime;
if ~isempty(ref)
    % plot current plan
    if check_if_plot_is_available(AH,'trajectory')
        AH.plot_data.trajectory.XData = Z(1,:) ;
        AH.plot_data.trajectory.YData = Z(2,:) ;
        AH.plot_data.trajectory.ZData = Z(3,:) ;
    else
        trajectory_data = plot3(Z(1,:),Z(2,:),Z(3,:),'k-',...
            'LineWidth',3) ;
        AH.plot_data.trajectory = trajectory_data ;
    end
end




end