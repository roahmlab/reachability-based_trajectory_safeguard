function [cost, dcost]=highway_cost( K,kv,kh,x_des,start_tic,timeout)
% find derivative with respect to variable, here it is k peak and ky
% k_pk = 9 ;
% kv= 10;
%get final pos, subtrack initial pose, get second norm.
[x,y,~] = predict_highway_position( K,kv,kh );
% dy = alpha/120*5*tpk.^4+beta/24*4*tpk.^3+gam/6*3*tpk.^2 +a0/2*2*tpk+v0;
%turns out 
cost = sum((x_des -[x;y]).^2);
dcost = zeros(2,1);
%dc/dx* dx/dpk
dcost(1) = -2*(x_des(1)-x)*1;
dcost(2) = -2*(x_des(2)-y)*19/48;
if toc(start_tic) > timeout
        error('Timed out while evaluating cost function!')
end
%Also need cost gradient with respect to K, not x and y

% del_v = 
% khh =kv*sin(kh);
% dely=ky-khh*tf;
% c31 =tf;
% delh = -khh;
% c121=1/c31^3*[-12 6*c31; 6*c31 -2*c31^2]*[dely;delh];
% c11= c121(1); c21= c121(2);
% f3= c11/6*t1.^3+c21/2*t1.^2+khh*t1 ;
%%

% accel = diff([0 f3]);
% f3=cumtrapz(T,accel.*F./max(F));
% hf = atan2(f3(end),f1(end));
% vy2= (f3(end)-f3(end-1))/(t1(end)-t1(end-1));
% dely2 = -f3(end);
% c32 = tf-tpk;
% delv2=-vy2;
% c122= 1/c32^3*[-12 6*c32; 6*c32 -2*c32^2]*[dely2;delv2];
% c12=c122(1);c22=c122(2);
% f4 = c12/6*t2.^3+c22/2*t2.^2+ vy2*t2+ f3(end);%ka*t2 ;
% f4= f3(end)*ones(1,101);
% f4_Z = f4;
% Z= [pos; f3_Z f4_Z(2:end)];
% 
% F= [F;f3 f4(2:end)];
% %%