function C = polygonsToCell(P)
% Given polygons in a 2-by-N array P, get each one into its own cell of the
% array C

% check that the first column and last column aren't nans
if isnan(P(1,1))
    P = P(:,2:end) ;
end

if isnan(P(1,end))
    P = P(:,1:end-1) ;
end

% get the number of polygons and the start indices of each one
N = size(P,2) ;
idx0 = 1:N ;
idxlog = isnan(P(1,:)) ;
idx0 = [1,(idx0(idxlog)+1)] ; % start indices
idx1 = [(idx0(2:end)-2),N] ; % end indices
NP = length(idx0) ;

% initialize the output
C = cell(1,NP) ;

% for each polygon, extract it and put it in C
for oidx = 1:NP
    i0 = idx0(oidx) ;
    i1 = idx1(oidx) ;
    p = P(:,i0:i1) ;
    C{oidx} = p ;
end

end