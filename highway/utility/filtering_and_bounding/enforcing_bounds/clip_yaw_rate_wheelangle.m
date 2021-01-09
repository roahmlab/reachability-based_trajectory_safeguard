function [w_filt,d_filt]=clip_yaw_rate_wheelangle(w,v,dmax,l)
%take in a vector of yaw rates and velocities,
%clip yaw rate and wheel angle based on kinematic assumption

d_filt=atan(l*w/v);

%reshape to columns
d_filt=d_filt(:);
w_filt=w(:);

L_v0=(v==0);

d_filt(L_v0)=0;
w_filt(L_v0)=0;

L_dmax=(abs(atan(l*w/v))>dmax);

d_filt(L_dmax)=sign(w)*dmax;
w_filt(L_dmax)=v/l.*tan(sign(w)*dmax);


   
end