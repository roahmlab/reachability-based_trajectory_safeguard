function wkstr = sub_z_and_t_into_w(wstr,ZT)
% wP = sub_z_into_w(wstr,ZT)
%
% INPUTS:
%   wstr  a structure representing the polynomial w, which must contain the
%         outputs from the msspoly function decomp, as the fields:
%           wvars   vector of free msspoly vars called k1,k2,z1,z2
%           wpows   m-by-4 array of powers of each term of w
%           wcoef   1-by-m array of the coefficients of w
%           zcols   1-by-2 vector of the columns of wpows for z1 and z2
%                   (which must be in the order [z1col z2col])
%           kcols   1-by-2 vector of the columns of wpows for k1 and k2 (in
%                   that order)
%           tcols   scalar of the column of wpows for t
%
%   ZT    a Nz-by-N vector of [Z;T]:
%             Z     Nz-by-N points in coordinates (z1,z2) to plug in to w
%                   and return an N-by-1 list of polynomials in k1 and k2
%             T     1-by-N points in t to plug in to w
%
% OUTPUTS:
%   wkstr a structure representing an N-by-1 list of polynomials in k1 and
%         k2, with fields wkpows, wkcoef, N, and kcols

Nz=length(wstr.zcols);

Z = ZT(1:Nz,:) ; T = ZT(Nz+1,:) ;

% extract info from w
wpows = wstr.wpows ;
wcoef = wstr.wcoef ;

% get the number of points to be evaluated
N = size(Z,2) ;

% get the powers of z1, z2, and t that will be used to evaluate w(t,z,.)
sub_pows = wpows(:,[wstr.zcols, wstr.tcols]) ;

% every point in P(1,:) needs to be raised to every power in sub_pows(:,1),
% and the same for P(2,:) wrt sub_pows(:,2) and P(3,:) wrt sub_pows(:,3)
P = sub_pows' ;

% get points in the format to the size: 3*Nterms x Npoints, where each pair
% of rows is alternating [z1;z2;t;z1;z2;t;...] to match the number of terms
% in the polynomial
ZTmat = repmat([Z ; T],size(P,2),1) ;

% get powers to the same size as the points, where each pair of rows
% corresponds to [z1;z2;t;...] and each column is repeated to be applied to
% each point (which is a column of [Z ; T])
Pmat = repmat(P(:),1,N) ;

% now raise the points to their corresponding powers
ZTmat = ZTmat.^Pmat ;

% collapse ZTmat by multiplying first rows against second rows against
% third rows, i.e. multiply the z1, z2, and t terms
if Nz == 2
    
    ZTmat = ZTmat(1:3:end,:).*ZTmat(2:3:end,:).*ZTmat(3:3:end,:);
    
else %if the size of z isn't 2 do something dumb
    
    for i=1:Nz+1
        if i==1
            ZTmattemp=ZTmat(1:Nz+1:end,:);
        else
            ZTmattemp=ZTmattemp.*(ZTmat(i:Nz+1:end,:));
        end
    end
    
    ZTmat = ZTmattemp;
    
end
% create a new coefficients matrix by replicating the current one Nterms
% times and multiplying each copy with the corresponding column in ZTmat
wkcoef = repmat(wcoef,N,1).*ZTmat' ;

% create the new powers matrix of the N polynomials in k1 and k2
wkpows = wpows(:,wstr.kcols) ;

%% COLLAPSING ROWS
% get all rows of wpows that are matched and sum the corresponding
% coefficients to create a row of coefficients that correspond to the
% reduced rows of powers

% sort wkpows first
[wkpowssorted,is] = sortrows(wkpows) ;
[wkpowsunique,~,ic] = unique(wkpowssorted,'rows','sorted') ;

% sort wkcoef's columns in the same order as wkpows' rows
wkcoef = wkcoef(:,is) ;

% sum the columns of wkcoef in groups based on how wkpows was made unique
% row by row (so identical rows of wkpows correspond to disparate columns
% of wkcoef, and these columns must be summed because they correspond to
% the terms in the polynomial with the same variables and powers)
summed_columns = splitapply(@(c) {sum(c,2)}, wkcoef, ic') ;

% concatenate each summed group of columns left to right
wkcoefcollapsed = cat(2,summed_columns{:}) ;

%% FINAL OUTPUT
% create the output structure
wkstr.wkcoef = wkcoefcollapsed ;
wkstr.wkpows = wkpowsunique ;
wkstr.N = N ;
wkstr.kcols = 1:length(wstr.kcols); % this is a formality
end

% function wkstr = sub_z_and_t_into_w(wstr,ZT)
% % wP = sub_z_into_w(wstr,ZT)
% %
% % INPUTS:
% %   wstr  a structure representing the polynomial w, which must contain the
% %         outputs from the msspoly function decomp, as the fields:
% %           wvars   vector of free msspoly vars called k1,k2,z1,z2
% %           wpows   m-by-4 array of powers of each term of w
% %           wcoef   1-by-m array of the coefficients of w
% %           zcols   1-by-2 vector of the columns of wpows for z1 and z2
% %                   (which must be in the order [z1col z2col])
% %           kcols   1-by-2 vector of the columns of wpows for k1 and k2 (in
% %                   that order)
% %           tcols   scalar of the column of wpows for t
% %
% %   ZT    a 3-by-N vector of [Z;T]:
% %             Z     2-by-N points in coordinates (z1,z2) to plug in to w
% %                   and return an N-by-1 list of polynomials in k1 and k2
% %             T     1-by-N points in t to plug in to w
% %
% % OUTPUTS:
% %   wkstr a structure representing an N-by-1 list of polynomials in k1 and
% %         k2, with fields wkpows, wkcoef, N, and kcols
% 
% Z = ZT(1:2,:) ; T = ZT(3,:) ;
% 
% % extract info from w
% wpows = wstr.wpows ;
% wcoef = wstr.wcoef ;
% 
% % get the number of points to be evaluated
% N = size(Z,2) ;
% 
% % get the powers of z1, z2, and t that will be used to evaluate w(t,z,.)
% sub_pows = wpows(:,[wstr.zcols, wstr.tcols]) ;
% 
% % every point in P(1,:) needs to be raised to every power in sub_pows(:,1),
% % and the same for P(2,:) wrt sub_pows(:,2) and P(3,:) wrt sub_pows(:,3)
% P = sub_pows' ;
% 
% % get points in the format to the size: 3*Nterms x Npoints, where each pair
% % of rows is alternating [z1;z2;t;z1;z2;t;...] to match the number of terms
% % in the polynomial
% ZTmat = repmat([Z ; T],size(P,2),1) ;
% 
% % get powers to the same size as the points, where each pair of rows
% % corresponds to [z1;z2;t;...] and each column is repeated to be applied to
% % each point (which is a column of [Z ; T])
% Pmat = repmat(P(:),1,N) ;
% 
% % now raise the points to their corresponding powers
% ZTmat = ZTmat.^Pmat ;
% 
% % collapse ZTmat by multiplying first rows against second rows against
% % third rows, i.e. multiply the z1, z2, and t terms
% ZTmat = ZTmat(1:3:end,:).*ZTmat(2:3:end,:).*ZTmat(3:3:end,:) ;
% 
% % create a new coefficients matrix by replicating the current one Nterms
% % times and multiplying each copy with the corresponding column in ZTmat
% wkcoef = repmat(wcoef,N,1).*ZTmat' ;
% 
% % create the new powers matrix of the N polynomials in k1 and k2
% wkpows = wpows(:,wstr.kcols) ;
% 
% %% COLLAPSING ROWS
% % get all rows of wpows that are matched and sum the corresponding
% % coefficients to create a row of coefficients that correspond to the
% % reduced rows of powers
% 
% % sort wkpows first
% [wkpowssorted,is] = sortrows(wkpows) ;
% [wkpowsunique,~,ic] = unique(wkpowssorted,'rows','sorted') ;
% 
% % sort wkcoef's columns in the same order as wkpows' rows
% wkcoef = wkcoef(:,is) ;
% 
% % sum the columns of wkcoef in groups based on how wkpows was made unique
% % row by row (so identical rows of wkpows correspond to disparate columns
% % of wkcoef, and these columns must be summed because they correspond to
% % the terms in the polynomial with the same variables and powers)
% summed_columns = splitapply(@(c) {sum(c,2)}, wkcoef, ic') ;
% 
% % concatenate each summed group of columns left to right
% wkcoefcollapsed = cat(2,summed_columns{:}) ;
% 
% %% FINAL OUTPUT
% % create the output structure
% wkstr.wkcoef = wkcoefcollapsed ;
% wkstr.wkpows = wkpowsunique ;
% wkstr.N = N ;
% wkstr.kcols = [1 2] ; % this is a formality
% end