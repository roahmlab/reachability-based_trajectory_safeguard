function [E,Ehi,Elo,E_all] = compute_2d_position_error(A,T_ref,Z_ref)
    % extract agent info
    T = A.time ;
    Z = A.state ;
    
    % get position trajectories
    Z = match_trajectories(T_ref,T,Z) ;
    h = Z(3,:);
    xy = reshape(Z(1:2,:),[],1);
    X = rotmat(h)*repmat(A.footprint_vertices,[length(h),1])+repmat(xy,[1,size(A.footprint_vertices,2)]);
%     X_plot=reshape(X,2,[]);
%     scatter(X_plot(1,:),X_plot(2,:));
%     A.plot();
    X_ref = repmat(reshape(Z_ref(1:2,:),[],1),[1,size(A.footprint_vertices,2)]) ;
    E_all = X- X_ref;
%     X = Z(1:2,:) ;
    Ehi = reshape(max(E_all,[],2),2,[]);
    Elo = reshape(min(E_all,[],2),2,[]);
    % compute position error
    E = reshape((Ehi+Elo)/2,2,[]) ;
end