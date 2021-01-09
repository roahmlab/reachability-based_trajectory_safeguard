%create a hexagonal grid in XY space, with points spaced evenly apart by
%point_spacing
%grid will at least cover more space than X_range and Y_range

function grid_points = hex_grid(X_range,Y_range,point_spacing,boundary)
if nargin<4
    boundary = 'none';
end

%find largest even number that is enough to span the distances
w=abs(point_spacing);
X_range = [min(X_range),max(X_range)];
Y_range = [min(Y_range),max(Y_range)];

if w>diff(X_range) && w>diff(Y_range)
    grid_points = [];
    return
end

X_range_expanded = 2/sqrt(3)*X_range+[-w,w];
Y_range_expanded = 2/sqrt(3)*Y_range+[-w,w];

num_w = max([ceil(diff(Y_range_expanded)/w)+(1-mod(ceil(diff(Y_range_expanded)/w),2)),...
    ceil(diff(X_range_expanded)/w)+(1-mod(ceil(diff(X_range_expanded)/w),2))]);

Xvec = X_range_expanded(1)+[0,1:num_w]*w;
Yvec = Y_range_expanded(1)+[0,1:num_w]*w;

[X,Y]=meshgrid(Xvec,Yvec);

n=size(Y,1);

X=sqrt(3)/2*X;
Y=Y+repmat([0,0.5],[n,n/2]);

X=reshape(X,[1,numel(X)]);
Y=reshape(Y,[1,numel(Y)]);

%get points that are in the X and Y range
if strcmp(boundary,'none')
    L=true(size(X));
elseif strcmp(boundary,'strict')
L=(X>X_range(1)&X<X_range(2)) & (Y>Y_range(1)&Y<Y_range(2));
else
    L=(X>=X_range(1)&X<=X_range(2)) & (Y>=Y_range(1)&Y<=Y_range(2));
end
grid_points =[X(L);Y(L)];

end