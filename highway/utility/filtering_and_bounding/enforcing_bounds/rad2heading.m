function d = rad2heading(h)
% Given an angle in radians h, convert it to a 'heading angle' d, in the
% interval (-pi,pi].
    d = mod(h+pi,2*pi) - pi ;
%     d = atan2(sin(h),cos(h)) ; % old version, slightly slower
end