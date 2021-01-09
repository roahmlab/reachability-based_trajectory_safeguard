function y_interp = interp_with_angles(x,y,q,angle_idxs,varargin)
%function expectes data to be in columns
%it uses the angle interpolations (radians) for indexs specified by
%angle_idxs
if numel(x)==1 && numel(q)==1 && all(q==x)

        y_interp = y;

else
    if nargin > 4
        y_interp = interp1(x(:),y,q(:),varargin{:});
        y_interp(:,angle_idxs)=interp_lon(x(:),y(:,angle_idxs),q(:),varargin{:});
    else
        y_interp = interp1(x(:),y,q(:));
        y_interp(:,angle_idxs)=interp_lon(x(:),y(:,angle_idxs),q(:));
    end
end

end