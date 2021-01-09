function [sym_list]=monomials_sym(sym_list_in,degree)
%input list of symbolic variables 1xn_sym, max total degree

%output a row vector off all monomial combinations of the input list from total degree 1 up to
%degree degree

%the function mss_asd is borrow from the Spotless toolbox
%https://github.com/spot-toolbox/spotless

nx=length(sym_list_in);
dd=mss_asd(nx,1:degree);
sym_list=prod(repmat(sym_list_in,[size(dd,1),1]).^dd,2)';
end
