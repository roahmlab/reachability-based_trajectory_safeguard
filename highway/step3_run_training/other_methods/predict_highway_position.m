function [x,y,v] = predict_highway_position( K,kv,kh)
ka = 0;
k_pk = K(1);
ky = K(2);
tpk = 2;
%%
delv1=k_pk-kv-ka*tpk;
c31 =tpk;
dela1 = -ka;
c121=1/c31^3*[-12 6*c31; 6*c31 -2*c31^2]*[delv1;dela1];
c11= c121(1); c21= c121(2);
v= c11/6*tpk.^3+c21/2*tpk.^2+ka*tpk + kv;
x = c11/24*tpk.^4+c21/6*tpk.^3+ka/2*tpk.^2 + kv*tpk;
%%
tlong = 4;
a0 = 0; 
v0= kv*sin(kh);
p0= 0;
del_p = ky-v0*tlong;
del_v = -v0 ;
consts =  1/tlong^5*[320 -120*tlong;
                -200*tlong 72*tlong^2;
                40*tlong^2 -12*tlong^3]*[del_p; del_v];
alpha = consts(1);
beta = consts(2);
gam = consts(3);

% t_vec = tlong;
y = alpha/120*tpk.^5+beta/24*tpk.^4+gam/6*tpk.^3 +a0/2*tpk.^2+v0*tpk+p0;
end