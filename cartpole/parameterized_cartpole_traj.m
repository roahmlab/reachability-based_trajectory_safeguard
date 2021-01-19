function [T, U, Z]=parameterized_cartpole_traj(k_pk,kv, ka)
%state 1: x state2:x dot state 3 theta, state 4 theta dot
if abs(k_pk) > 5
    k_pk = 5*sign(k_pk);
end
tf= 0.3;
tpk = 0.1;
t1 =linspace(0, tpk,101);
t2 =linspace(0, tf-tpk,101);

delv1=k_pk-kv-ka*tpk;
c31 =tpk;
dela1 = -ka;

% compute spline parameters
% [ax, bx, cx] = single_axis_params(Dv,Da,t_peak);
% 
% % write velocity
% dp = (ax/24).*t.^4 + (bx/6).*t.^3 + (cx/2).*t.^2 + (a_0).*t + v_0;
c121=1/c31^3*[-12 6*c31; 6*c31 -2*c31^2]*[delv1;dela1];
c11= c121(1); c21= c121(2);
f1= c11/6*t1.^3+c21/2*t1.^2+ka*t1 + kv;

delv2 = -k_pk;
c32 = tf-tpk;
 dela2=0;
c122= 1/c32^3*[-12 6*c32; 6*c32 -2*c32^2]*[delv2;dela2];
c12=c122(1);c22=c122(2);
f2 = c12/6*t2.^3+c22/2*t2.^2+ k_pk;%ka*t2 ;

T = [t1 t2(2:end)+tpk];
U = [f1 f2(2:end)];
Z = [cumsum(f1)*(t1(2)-t1(1))];
Z = [Z Z(end)+cumsum(f2(2:end))*(t2(2)-t2(1))];
end
