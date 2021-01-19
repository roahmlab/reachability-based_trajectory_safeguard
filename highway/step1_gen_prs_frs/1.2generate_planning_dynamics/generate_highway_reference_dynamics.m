function [] = generate_highway_reference_dynamics()
% computing and saving the planning model for FRS
% computation. see parameterized_traj_1d_with_brake.m for more details
%t1  -------t2 --------t3------ t
%f1 ------f2(f21) ---(f22)----- v
%f3(f31)--f3(f32)-----f4------- dy/dt
% write overall dynamics
%FRS states:
%delta   0    des   ino  ini  des  
%v       y    k_pk, kv,  kh,  ky,  t
% note to reduce number of uncertainty, v only account for delta velocity, so initial condition uncertainty only have to appear in kv
tpk = 2; 
tf = 6;
ka = 0;
%% initial condition params
syms v y k_pk kv kh ky t1 
%% Compute "to peak" dynamics
delv1=k_pk-kv-ka*tpk;
c31 =tpk;
dela1 = -ka;
c121=1/c31^3*[-12 6*c31; 6*c31 -2*c31^2]*[delv1;dela1];
c11= c121(1); c21= c121(2);
f1= c11/6*t1.^3+c21/2*t1.^2+ka*t1 + kv;

tmid=4;
del_v = -kv*sin(kh) ;
del_p = ky + del_v*tmid;

consts =  1/tmid^5*[320 -120*tmid(end);
                -200*tmid(end) 72*tmid(end)^2;
                40*tmid(end)^2 -12*tmid(end)^3]*[del_p; del_v];
alpha = consts(1);
beta = consts(2);
gam = consts(3);

t_vec = t1;
%derivative of the y_reference, integrate to get command U or F, here it needs
%vy to propogate forward properly in CORA
f31 = alpha/24*t_vec.^4+beta/6*t_vec.^3+gam/2*t_vec.^2 -del_v;


% x = [v; v_0; a_0; v_peak; t];
x = [v; y; k_pk; kv; kh;ky; t1];

 syms tdummy udummy
% 
%  dx = [dp; 0; 0; 0; 1];
dx = [f1;f31; 0; 0;0;0; 1 ]
matlabFunction(dx, 'File', 'highway_toPeak', 'vars', {tdummy x udummy});

%% Compute "to mid" dynamics
syms t2

t2_plus2 = t2+tpk*ones(size(t2));
f32 = alpha/24*t2_plus2.^4+beta/6*t2_plus2.^3+gam/2*t2_plus2.^2 +a0*t2_plus2+v0;

delv2 = -k_pk;
c32 = 4;
 dela2=0;
c122= 1/c32^3*[-12 6*c32; 6*c32 -2*c32^2]*[delv2;dela2];
c12=c122(1);c22=c122(2);
f21 = c12/6*t2.^3+c22/2*t2.^2+ k_pk;

x = [v; y; k_pk; kv; kh;ky; t2];
dx = [f21;f32; 0; 0;0;0; 1 ]
matlabFunction(dx, 'File', 'highway_mid', 'vars', {tdummy x udummy});

%% Compute "to stop" dynamics
syms t3

delv2 = -k_pk;
c32 = 4;
 dela2=0;
c122= 1/c32^3*[-12 6*c32; 6*c32 -2*c32^2]*[delv2;dela2];
c12=c122(1);c22=c122(2);
t3_plus2 = t3+tpk*ones(size(t3));
f22 = c12/6*t3_plus2.^3+c22/2*t3_plus2.^2+ k_pk;

x = [v; y; k_pk; kv; kh;ky; t3];
dx = [f22;0; 0; 0;0;0; 1 ];
matlabFunction(dx, 'File', 'highway_toStop', 'vars', {tdummy x udummy});

end

   
