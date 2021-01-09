function dd=mss_asd(m,p)
%this function is from the spotless toolbox
%https://github.com/spot-toolbox/spotless

% function dd=mss_asd(m,p)
%
% INPUTS:
%   m   -  a positive integer 
% for a positive integer m and a vector p of different integers
% rows of dd (N-by-m) represent all ways of representing elements of p
% as sums of m non-negative integers

if nargin<2, error('2 inputs required'); end
if ~isa(m,'double'), error('1st input not a double'); end
m=max(1,round(m(1)));
if ~isa(p,'double'), error('2nd input not a double'); end
p=p(:);
if ~all(p==max(0,round(p))), error('2nd input not integer>=0'); end
p=sort(p);
np=length(p);
if any([-1;p(1:np-1)]==p), error('2nd input has repeated entries'); end
smx=max(p);
if smx==0, dd=zeros(1,m); return; end
N=0;
for i=1:np, N=N+nchoosek(p(i)+m-1,m-1); end
dd=zeros(N,m);
d=[zeros(1,m-1) p(1)];       % minimal solution
t=1;                         % current row in dd
i=m;                         % last non-zero element in d
k=1;                         % k such that s=p(k)
while 1,
    dd(t,:)=d;           % write the result
    t=t+1;               % next dd row
    if t>N, break; end
    if k==np,               % reached the maximum  
        s=p(k)-d(i)+1;       % next sum
        k=1+length(p(p<s));  % next k
        d(i)=0;              % zero out current position
        d(i-1)=d(i-1)+1;     % up the previous one
        d(m)=p(k)-s;
        if d(m)>0, 
            i=m; 
        else
            i=i-1;
        end
    else
        d(m)=d(m)+p(k+1)-p(k);
        k=k+1;
        i=m;
    end
end
