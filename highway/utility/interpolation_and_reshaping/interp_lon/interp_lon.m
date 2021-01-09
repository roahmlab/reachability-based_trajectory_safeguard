function out = interp_lon(x,lon,xq,varargin)

% INTERP_LON interpolates a set of longitude angles (in deg)
%
% Usage: out = interp_lon(x,lon,xq)
%
% x and lon are vectors of length N.  function evalutes longitude 
% (in deg -180..180) at points xq using unwrap and interp1
%
% to specify interpolation method used in interp1, use
% out = interp_lon(x,lon,xq,METHOD)
%
% Written by D.G. Long, 27 Nov 2017 

%modified by Sean on 6 Fed 2019 to work with radians

ulon=unwrap(lon);
if nargin>3
  out=interp1(x,ulon,xq,varargin{1});
else
  out=interp1(x,ulon,xq);
end
out=mod(out,2*pi);
out(out>pi)=out(out>pi)-2*pi;

end