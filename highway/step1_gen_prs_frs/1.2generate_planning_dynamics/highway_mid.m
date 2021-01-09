function dx = highway_mid(tdummy,in2,udummy)
%HIGHWAY_MID
%    DX = HIGHWAY_MID(TDUMMY,IN2,UDUMMY)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    27-Jul-2020 12:15:03

k_pk = in2(3,:);
kh = in2(5,:);
kv = in2(4,:);
ky = in2(6,:);
t2 = in2(7,:);
t3 = sin(kh);
t4 = t2+2.0;
dx = [k_pk-k_pk.*t2.^2.*(3.0./1.6e+1)+(k_pk.*t2.^3)./3.2e+1;kv.*t3+t4.^2.*(ky.*(5.0./1.6e+1)-kv.*t3.*(7.0./8.0))-t4.^3.*(ky.*(2.5e+1./1.92e+2)-(kv.*t3)./3.0)+t4.^4.*(ky.*(5.0./3.84e+2)-kv.*t3.*(2.5e+1./7.68e+2));0.0;0.0;0.0;0.0;1.0];
