function plot_w_contour_filled(wz,z,z0,x0,y0,Dx,Dy,face_color,edge_color,face_alpha,N,R)
    % plot FRS
    if ~exist('N','var')
        N = 250 ;
    end
    if ~exist('R','var')
        R = 0.8;
    end

    xvec = linspace(-R,R,N) ;
    yvec = linspace(-R,R,N) ;
    [X,Y] = meshgrid(xvec, yvec) ;
    ZZ = [X(:) Y(:)]' ;
    F = reshape(full(msubs(wz,z(1:2), ZZ)),N,N) ;

    [hh]=contourc(xvec,yvec,F,[1 1]);             

    idxs=find((hh(1,:)==1)&(mod(hh(2,:),1)==0));
    x=[];
    y=[];
    
    for idx=idxs
        x=[x,hh(1,idx+1:idx+hh(2,idx))];
        y=[y,hh(2,idx+1:idx+hh(2,idx))];
    end
        xshift=(Dx*x-x0)*cos(z0(3))-sin(z0(3))*(Dy*y-y0);
        yshift=(Dy*y-y0)*cos(z0(3))+sin(z0(3))*(Dx*x-x0);

        FRS_contour=z0(1:2)+[xshift;yshift];
        hold on
%         plot(FRS_contour(1,:),FRS_contour(2,:),'Color',color,'LineWidth',1.5,'LineStyle',linestyle)
        patch(FRS_contour(1,:)',FRS_contour(2,:)',face_color,...
             'EdgeColor',edge_color,'FaceAlpha',face_alpha)
    
end