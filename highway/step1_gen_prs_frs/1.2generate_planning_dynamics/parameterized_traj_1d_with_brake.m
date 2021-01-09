% syms t k_pk kv tf tpk
%kv - 1 < k_pk< kv+ 0.5
function [T,F,Z]=parameterized_traj_1d_with_brake(k_pk, kv,kh,ky)
% k_pk = 9 ;
% kv= 10;
ka = 0;
% 
tf= 6;
tpk = 2;
t1 =linspace(0, tpk,101);
t2 =linspace(0, tf-tpk,201);
%%
delv1=k_pk-kv-ka*tpk;
c31 =tpk;
dela1 = -ka;
c121=1/c31^3*[-12 6*c31; 6*c31 -2*c31^2]*[delv1;dela1];
c11= c121(1); c21= c121(2);
f1= c11/6*t1.^3+c21/2*t1.^2+ka*t1 + kv;
% f1_Z = c11/24*t1.^4+c21/6*t1.^3+ka/2*t1.^2 + kv*t1;
%%
delv2 = -k_pk;
c32 = tf-tpk;
 dela2=0;
c122= 1/c32^3*[-12 6*c32; 6*c32 -2*c32^2]*[delv2;dela2];
c12=c122(1);c22=c122(2);

% c12 = -3*delv2; c22= 3*c32*delv2;

f2 = c12/6*t2.^3+c22/2*t2.^2+ k_pk;%ka*t2 ;
% f2_Z = c12/24*t2.^4+c22/6*t2.^3+ka/2*t2.^2 + kv*t2+f1_Z(end)*ones(size(t2));

T= [t1 t2(2:end)+tpk];
F= [f1 f2(2:end)];

pos = cumsum(F)*(t1(2)-t1(1));
%%
tlong=T(1:201);
a0 = 0; 
v0= kv*sin(kh);
p0= 0;
del_p = ky-v0*tlong(end);
del_v = -v0 ;
consts =  1/tlong(end)^5*[320 -120*tlong(end);
                -200*tlong(end) 72*tlong(end)^2;
                40*tlong(end)^2 -12*tlong(end)^3]*[del_p; del_v];
alpha = consts(1);
beta = consts(2);
gam = consts(3);

t_vec = tlong;
f3 = alpha/120*t_vec.^5+beta/24*t_vec.^4+gam/6*t_vec.^3 +a0/2*t_vec.^2+v0*t_vec+p0;
f3_Z = f3;
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
f4= f3(end)*ones(1,101);
f4_Z = f4;
Z= [pos; f3_Z f4_Z(2:end)];

F= [F;f3 f4(2:end)];
% %%
% figure(11);clf;hold on;
% plot(T,F(1,:))
% % plot(t2+tpk,f2)
% figure(12);clf;hold on;
% plot(T(1:end-1),diff(F(1,:))./diff(T));
% 
% figure(13);clf;hold on;
% plot(T,F(2,:));


%this thing here cannot deal with dynamic where the first 0.5s the desired
%velocity needs to be something else, it has to be piecewise, tanh a better
%option.
%vx + (2*vx*t^3)/(tf - tpk)^3 - (vx*t^2*(6*tf - 6*tpk))/(2*(tf - tpk)^3)