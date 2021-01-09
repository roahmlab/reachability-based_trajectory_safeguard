function handle=plot_w_contour(wz,z,z0,x0,y0,Dx,Dy,N,range,varargin)
    % plot FRS
    
    if nargin<10
        varargin = [];
        
        if nargin<9
            range = 1;
            if nargin<8
                N = 250;
                 
            end
        end
    end
   if isempty(range)
       range=1;
   end
   if isempty(N)
       N=250;
   end

    if numel(range)==1
        xvec = linspace(-range,range,N) ;
        yvec = linspace(-range,range,N) ;
    else
        xvec = linspace(range(1,1),range(1,2),N) ;
        yvec = linspace(range(2,1),range(2,2),N) ;
    end
    
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
        hold on
        
        if ~isempty(varargin)
            
            if numel(varargin)==1
                if iscell(varargin{1})
                     varargin = varargin{:};
                end
            end
            
            if nargout==1
                handle=plot(FRS_contour(1,:),FRS_contour(2,:),varargin{:});
            else
                plot(FRS_contour(1,:),FRS_contour(2,:),varargin{:})
            end
        else
            if nargout==1
                handle=plot(FRS_contour(1,:),FRS_contour(2,:));
            else
                plot(FRS_contour(1,:),FRS_contour(2,:))
            end
            
        end
    end
    % colormap(color)
    % contour(xvec,yvec,F,[1 1])
end