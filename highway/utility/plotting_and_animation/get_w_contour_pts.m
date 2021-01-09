function FRS_contour = get_w_contour_pts(wz,z,z0,x0,y0,Dx,Dy,space_bound,N)

    if nargin <9
        N=100;
    if nargin <8
        space_bound = 0.9;
    end
    end
    

    xvec = linspace(-space_bound,space_bound,N) ;
    yvec = linspace(-space_bound,space_bound,N) ;
    [X,Y] = meshgrid(xvec, yvec) ;
    ZZ = [X(:) Y(:)]' ;
    F = reshape(full(msubs(wz,z(1:2), ZZ)),N,N) ;

    [hh]=contourc(xvec,yvec,F,[1 1]);             

    idxs=find((hh(1,:)==1)&(mod(hh(2,:),1)==0));

    for idx=idxs
        x=hh(1,idx+1:idx+hh(2,idx));
        y=hh(2,idx+1:idx+hh(2,idx));

        xshift=(Dx*x-x0)*cos(z0(3))-sin(z0(3))*(Dy*y-y0);
        yshift=(Dy*y-y0)*cos(z0(3))+sin(z0(3))*(Dx*x-x0);

        FRS_contour=z0(1:2)+[xshift;yshift];

    end
end