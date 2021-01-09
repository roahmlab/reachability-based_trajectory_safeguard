function [v] = get_derivative(x,t)
N=length(t);
v=NaN(N,1);
v(1)=(x(2)-x(1))/(t(2)-t(1));
for i=2:N-1
    v(i)=(x(i+1)-x(i))/(t(i+1)-t(i));
end
v(N)=(x(N)-x(N-1))/(t(N)-t(N-1));

end

