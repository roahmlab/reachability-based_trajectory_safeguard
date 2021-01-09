% syms t k_pk kv tf tpk
%kv - 1 < k_pk< kv+ 0.5
function [F]=parameterized_traj_1d(k_pk, kv,ka, tf)
% k_pk = 9 ;
% kv= 10;
% ka = 0;

% tf= 5;
 tpk = tf;
t1 =linspace(0, tpk,200);
% t2 =linspace(0, tf-tpk);
%%
delv1=k_pk-kv-ka*tpk;
c31 =tpk;
dela1 = -ka;
c121=1/c31^3*[-12 6*c31; 6*c31 -2*c31^2]*[delv1;dela1];
c11= c121(1); c21= c121(2);
f1= c11/6*t1.^3+c21/2*t1.^2+ka*t1 + kv;
%%
% delv2 = -k_pk;
% c32 = tf-tpk;
% dela2=0;
% c122= 1/c32^3*[-12 6*c32; 6*c32 -2*c32^2]*[delv2;dela2];
% c12=c122(1);c22=c122(2);
% 
% % c12 = -3*delv2; c22= 3*c32*delv2;
% 
% f2 = c12/6*t2.^3+c22/2*t2.^2+ k_pk;%ka*t2 ;
T= t1;
F= f1;
% % %%
% figure(1);clf;hold on;
% plot(t1,f1)
% plot(t2+tpk,f2)
% figure(2);clf;hold on;
% plot(t1(1:end-1),diff(f1)./diff(t1));
% t2_real = t2+tpk;
% plot(t2_real(1:end-1),diff(f2)./diff(t2_real))

%this thing here cannot deal with dynamic where the first 0.5s the desired
%velocity needs to be something else, it has to be piecewise, tanh a better
%option.
%vx + (2*vx*t^3)/(tf - tpk)^3 - (vx*t^2*(6*tf - 6*tpk))/(2*(tf - tpk)^3)