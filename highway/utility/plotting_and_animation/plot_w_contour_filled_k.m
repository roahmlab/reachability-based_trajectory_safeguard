function plot_w_contour_filled_k(wk,k,k_range,face_color,edge_color,face_alpha)
    % plot FRS
    N = 100 ;
   
        

    k1vec = linspace(-1,1,N) ;
    k2vec = linspace(-1,1,N) ;
    [K1,K2] = meshgrid(k1vec, k2vec) ;
    KK = [K1(:) K2(:)]' ;
    
    F = reshape(full(msubs(wk,k(1:2), KK)),N,N) ;
    
    cmap = colormap;
    map = face_color;
    colormap(map)
    
      contourf(diff(k_range(1,:))/2*(k1vec+1)+k_range(1,1),diff(k_range(2,:))/2*(k2vec+1)+k_range(2,1),F,[1 1],'Color',edge_color);
     
     
% 
%    [hh]=contourc(k1vec,k2vec,F,linspace(1,max(max(F)),10));             
% 
%     idxs=find((hh(1,:)==1)&(mod(hh(2,:),1)==0));
%     start=true;
%     for idx=idxs
%         if idx<size(hh,2)
%         x=hh(1,idx+1:idx+hh(2,idx));
%         y=hh(2,idx+1:idx+hh(2,idx));
%         
%         x=diff(k_range(1,:))/2*(x+1)+k_range(1,1);
%         
%         y=diff(k_range(2,:))/2*(y+1)+k_range(2,1);
%         
%         FRS_contour=[x;y];
%         hold on
%         
%         if start
%         
%         patch(FRS_contour(1,:),FRS_contour(2,:),face_color,...
%              'EdgeColor',edge_color,'LineWidth',1.5,'FaceAlpha',face_alpha)
%         else
%                 patch(FRS_contour(1,:),FRS_contour(2,:),face_color,...
%              'EdgeColor','none','LineWidth',1.5,'FaceAlpha',face_alpha)
%         end
%         end
%     end
end