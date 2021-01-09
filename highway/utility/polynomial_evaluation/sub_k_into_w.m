function wstr = sub_k_into_w(wstr,k,kcolidx)

% extract info from w
wpows = wstr.wpows ;
wcoef = wstr.wcoef ;
kcols = wstr.kcols ;
zcols = wstr.zcols ;

% raise the point x to its corresponding power in each monomial
col = kcols(kcolidx) ;
sub_pows = wpows(:,col) ;
Xmat = k.^(sub_pows') ;

% multiply the coefficients by the resulting evaluation
wcoef = wcoef.*Xmat ;

% create the new powers matrix, and adjust the column indices
wpows(:,col) = [] ;
kcols(kcolidx) = [] ;
kcols(kcols > col) = kcols(kcols > col) - 1 ;
zcols(zcols > col) = zcols(zcols > col) - 1 ;


%% COLLAPSING ROWS
% % sort wkpows first
% [wpowssorted,is] = sortrows(wpows) ;
% [wpowsunique,~,ic] = unique(wpowssorted,'rows','sorted') ;
% 
% % sort wkcoef's columns in the same order as wkpows' rows
% wcoef = wcoef(:,is) ;
% 
% % sum the columns of wkcoef in groups based on how wkpows was made unique
% % row by row (so identical rows of wkpows correspond to disparate columns
% % of wkcoef, and these columns must be summed because they correspond to
% % the terms in the polynomial with the same variables and powers)
% summed_columns = splitapply(@(c) {sum(c,2)}, wcoef, ic') ;
% 
% % concatenate each summed group of columns left to right
% wcoefcollapsed = cat(2,summed_columns{:}) ;

%% FINAL OUTPUT
% create the output structure
% wstr.wcoef = wcoefcollapsed ;
% wstr.wpows = wpowsunique ;
wstr.wcoef = wcoef ;
wstr.wpows = wpows ;
wstr.kcols = kcols ;
wstr.zcols = zcols ;
end