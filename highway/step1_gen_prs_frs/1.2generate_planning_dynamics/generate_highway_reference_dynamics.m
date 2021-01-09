function [] = generate_highway_reference_dynamics()
% computing and saving the quadcopter dynamics functions for FRS
% computation.
% v_peak is the only control parameter
% compute two dynamics functions: one for the time until v_peak, and one
% from v_peak to stopping behavior.

%FRS states:
%v y k_pk, kv,kh,ky, t

tpk = 2; 
tf = 6;
% t_to_stop = t_total - t_peak; % time from peak to stop
% p0 = 0; % assume initial position is set as the origin.
% a_peak = 0; % 0 acceleration at peak speed
% v_f = 0; % 0 final speed 
% a_f = 0; % 0 final acceleration.
%% User defined params
% k_pk = 0;
% ky = 0;
 ka = 0;
%% initial condition params
syms v y k_pk kv kh ky t1 
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
% write peak position
% p_plan = (ax/120).*t_plan.^5 + (bx/24).*t_plan.^4 + (cx/6).*t_plan.^3 + (a_0(1)/2).*t_plan.^2 + v_0(1).*t_plan;
% p_peak = (ax/120).*t_peak.^5 + (bx/24).*t_peak.^4 + (cx/6).*t_peak.^3 + (a_0(1)/2).*t_peak.^2 + v_0(1).*t_peak;
%%
tmid=4;
% v0= kv*sin(kh);
del_v = -kv*sin(kh) ;
del_p = ky + del_v*tmid;

consts =  1/tmid^5*[320 -120*tmid(end);
                -200*tmid(end) 72*tmid(end)^2;
                40*tmid(end)^2 -12*tmid(end)^3]*[del_p; del_v];
alpha = consts(1);
beta = consts(2);
gam = consts(3);

t_vec = t1;
%derivative of the y_reference, integrate to get command U, here it needs
%vy to propogate forward
f31 = alpha/24*t_vec.^4+beta/6*t_vec.^3+gam/2*t_vec.^2 -del_v;
%t1  -------t2 --------t3------
%f1 ------f2(f21) ---(f22)-----
%f3(f31)--f3(f32)-----f4-------
% write overall dynamics

% x = [v; v_0; a_0; v_peak; t];
x = [v; y; k_pk; kv; kh;ky; t1];

 syms tdummy udummy
% 
%  dx = [dp; 0; 0; 0; 1];
dx = [f1;f31; 0; 0;0;0; 1 ]
matlabFunction(dx, 'File', 'highway_toPeak', 'vars', {tdummy x udummy});

%% Compute "to stop" dynamics
syms t2

t2_plus2 = t2+tpk*ones(size(t2));
f32 = alpha/24*t2_plus2.^4+beta/6*t2_plus2.^3+gam/2*t2_plus2.^2 +a0*t2_plus2+v0;

delv2 = -k_pk;
c32 = 4;
 dela2=0;
c122= 1/c32^3*[-12 6*c32; 6*c32 -2*c32^2]*[delv2;dela2];
c12=c122(1);c22=c122(2);
f21 = c12/6*t2.^3+c22/2*t2.^2+ k_pk;%ka*t2 ;

x = [v; y; k_pk; kv; kh;ky; t2];
% 
%  dx = [dp; 0; 0; 0; 1];
dx = [f21;f32; 0; 0;0;0; 1 ]
matlabFunction(dx, 'File', 'highway_mid', 'vars', {tdummy x udummy});
%%
syms t3

delv2 = -k_pk;
c32 = 4;
 dela2=0;
c122= 1/c32^3*[-12 6*c32; 6*c32 -2*c32^2]*[delv2;dela2];
c12=c122(1);c22=c122(2);
t3_plus2 = t3+tpk*ones(size(t3));
f22 = c12/6*t3_plus2.^3+c22/2*t3_plus2.^2+ k_pk;%ka*t2 ;

x = [v; y; k_pk; kv; kh;ky; t3];
dx = [f22;0; 0; 0;0;0; 1 ];
matlabFunction(dx, 'File', 'highway_toStop', 'vars', {tdummy x udummy});


%UNDER CONSTRUCTION
% for each axis, compute the change in velocity/accel
% Dv = v_f - v_peak - a_peak*t_to_stop ;
% Da = a_f - a_peak ;
% 
% % compute spline parameters
% [ax, bx, cx] = single_axis_params(Dv,Da,t_to_stop);
% 
% % write velocity
% dp = (ax/24).*t.^4 + (bx/6).*t.^3 + (cx/2).*t.^2 + (a_peak).*t + v_peak;
% 
% % write overall dynamics
% x = [p; v_0; a_0; v_peak; t];
% syms tdummy udummy
% 
% dx = [dp; 0; 0; 0; 1];
% matlabFunction(dx, 'File', 'dyn_quadrotor_toStop', 'vars', {tdummy x udummy});
% 
% % compute the final position p_f as a function of v_peak
% p_f = (ax/120).*t_to_stop.^5 + (bx/24).*t_to_stop.^4 + (cx/6).*t_to_stop.^3 + (a_peak/2).*t_to_stop.^2 + v_peak.*t_to_stop + p_peak;
% 
% % save position functions:
% matlabFunction(p_plan, 'File', 'pos_quadrotor_plan', 'vars', {v_0 a_0 v_peak});
% matlabFunction(p_peak, 'File', 'pos_quadrotor_peak', 'vars', {v_0 a_0 v_peak});
% matlabFunction(p_f, 'File', 'pos_quadrotor_final', 'vars', {v_0 a_0 v_peak});

end

   
