function [v,eul,R] = getlocal_velocites(x,y,z,q,t)
N=length(t);
v=NaN(N,3);
R=quat2rotm(q);
eul=quat2eul(q,'XYZ');

dxdt=get_5point_derivative(x,t);
dydt=get_5point_derivative(y,t);
dzdt=get_5point_derivative(z,t);

for i=1:N
    v(i,:)=R(:,:,i)'*[dxdt;dydt;dzdt];
end

end

