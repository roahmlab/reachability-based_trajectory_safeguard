% Last Edited: Jan 2021
% Yifei Simon Shao
% This code generate the reference trajectory for car example in paper
% takes in desired velocity k_pk, initial velocity kv, initial heading kh
% and desired latitudinal offset ky.
% output is T vector for time [0,6]s, 2s to peak velocity, 4s to desired y, 6 seconds to stop. 
% F is desired control input in each dimension for the car controller
% Z is the reference state that may or may not be used by the controller
% syms t k_pk kv tf tpk
% this reference trajectory was proposed here for drones and adapted for our purpose:
% https://ieeexplore.ieee.org/document/7299672
function [T,F,Z]=parameterized_traj_1d_with_brake(k_pk, kv,kh,ky)
% k_pk = 9 ;
% kv= 10;
ka = 0; %assuming 0 initial acceleration, this will make the motion more jerky but does reduce the number of parameters to 2 
tf= 6; % 2 + 2 + 2
tpk = 2;
t1 =linspace(0, tpk,101);
t2 =linspace(0, tf-tpk,201);
%% reference velocity to desired velocity
delv1=k_pk-kv-ka*tpk;
c31 =tpk;
dela1 = -ka;
c121=1/c31^3*[-12 6*c31; 6*c31 -2*c31^2]*[delv1;dela1];
c11= c121(1); c21= c121(2);
f1= c11/6*t1.^3+c21/2*t1.^2+ka*t1 + kv; 
%% reference velocity to 0
delv2 = -k_pk; 
c32 = tf-tpk;
 dela2=0;
c122= 1/c32^3*[-12 6*c32; 6*c32 -2*c32^2]*[delv2;dela2];
c12=c122(1);c22=c122(2);
f2 = c12/6*t2.^3+c22/2*t2.^2+ k_pk;

T= [t1 t2(2:end)+tpk];
F= [f1 f2(2:end)];

pos = cumsum(F)*(t1(2)-t1(1)); % integrate velocity to get position
%% lat reference, this controller takes in desired lateral position as input so F will be position
tlong=T(1:201); % 0-4 seconds to desired y position
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

f4= f3(end)*ones(1,101); % maintain final y
f4_Z = f4;
Z= [pos; f3_Z f4_Z(2:end)]; % position in y is same as control input F

F= [F;f3 f4(2:end)];
