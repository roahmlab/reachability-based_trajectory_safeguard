function f_out = dynamic_fun_world_to_scaled_box(f,t,x,x_scale,x_offset,T)
if nargin<6
    T = 1;
end

N = size(f,2);

f_out = msspoly(zeros(size(f)));

for i = 1:N
f_out(:,i) = (T./x_scale) .* subs(f(:,i),[t;x],[T*t;x_scale.*x-x_offset]);
end

end
