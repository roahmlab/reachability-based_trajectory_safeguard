function plot_wt_contours_filled(wtz,t_mss,z_mss,t_vec,z0,x0,y0,Dx,Dy,start_color,end_color,edge_color,fig_num)
    % parse inputs
    if nargin < 13
        fig_num = 1 ;
    if nargin < 12
        edge_color = 'none' ;
    if nargin < 11
        end_color = [1 0.3 0.3] ;
    if nargin < 10
        start_color = [1 0.9 0.9] ;
    end
    end
    end
    end
    
    N = 100 ;
    
    for t = t_vec
    
    wz = msubs(wtz,t_mss,t);
    
  
        
    xvec = linspace(-0.9,0.9,N) ;
    yvec = linspace(-0.9,0.9,N) ;
    [X,Y] = meshgrid(xvec, yvec) ;
    ZZ = [X(:) Y(:)]' ;
    F = reshape(full(msubs(wz,z_mss(1:2), ZZ)),N,N) ;

    [hh]=contourc(xvec,yvec,F,[1 1]);             

    idxs=find((hh(1,:)==1)&(mod(hh(2,:),1)==0));
    if ~isempty(idxs)
    for idx=idxs
        x=hh(1,idx+1:idx+hh(2,idx));
        y=hh(2,idx+1:idx+hh(2,idx));

        xshift=(Dx*x-x0)*cos(z0(3))-sin(z0(3))*(Dy*y-y0);
        yshift=(Dy*y-y0)*cos(z0(3))+sin(z0(3))*(Dx*x-x0);

        FRS_contour=z0(1:2)+[xshift;yshift];
    end
    
       color_cur = interp1([min(t_vec);max(t_vec)],[start_color;end_color],t) ;
       
       figure(fig_num)
       
       patch(FRS_contour(1,:)',FRS_contour(2,:)',color_cur,'EdgeColor',edge_color)
    end
end
end