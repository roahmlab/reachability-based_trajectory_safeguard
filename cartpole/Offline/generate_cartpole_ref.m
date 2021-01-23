function [] = generate_cartpole_ref()
% computing and saving the quadcopter dynamics functions for FRS
tpk = 0.1; 
tf = 0.3;
% t_to_stop = t_total - t_peak; % time from peak to stop
% p0 = 0; % assume initial position is set as the origin.
% a_peak = 0; % 0 acceleration at peak speed
% v_f = 0; % 0 final speed 
% a_f = 0; % 0 final acceleration.
%% User defined params
% k_pk = 0;
%% initial condition params
syms v k_pk kv ka t1
%% Compute "to peak" dynamics

% syms p v_0 a_0 v_peak t; % variables for FRS: position, time. parameters for FRS: initial speed, initial acceleration, peak speed

% compute change in velocity/accel for each axis
% Dv = v_peak - v_0 - a_0*t_peak ;
% Da = a_peak - a_0 ;
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

x = [v; k_pk; kv;ka; t1];

 syms tdummy udummy
% 
%  dx = [dp; 0; 0; 0; 1];
dx = [f1;0;0;0; 1 ];

matlabFunction(dx, 'File', 'cartpole_toPeak', 'vars', {tdummy x udummy});

syms t2;

delv2 = -k_pk;
c32 = tf-tpk;
dela2=0;
c122= 1/c32^3*[-12 6*c32; 6*c32 -2*c32^2]*[delv2;dela2];
c12=c122(1);c22=c122(2);
f2 = c12/6*t2.^3+c22/2*t2.^2+ k_pk;%ka*t2 ;

x = [v; k_pk; kv;ka; t2];

% 
%  dx = [dp; 0; 0; 0; 1];
dx = [f2;0;0;0; 1 ];

matlabFunction(dx, 'File', 'cartpole_toStop', 'vars', {tdummy x udummy});
end
