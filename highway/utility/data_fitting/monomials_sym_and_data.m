function [sym_list,data]=monomials_sym_and_data(sym_list_in,data_in,degree)
%input: sym_list_in: list of symbolic variables 1xn_sym, max total degree
%        data_in: list of data for symbolic variables n_dataxn_terms, max total degree
%output: sym_list a row vector off all monomial combinations of the input list from total degree 1 up to
%          data: n_data row vectors off all monomial combinations of the input list from total degree 1 up to

%the function mss_asd is borrow from the Spotless toolbox
%https://github.com/spot-toolbox/spotless

nx=length(sym_list_in);
dd=mss_asd(nx,1:degree);
sym_list=prod(repmat(sym_list_in,[size(dd,1),1]).^dd,2)';

%output data length(sym_list) x ndata
sz=size(data_in);
if sz(1)==length(sym_list_in)
     data_in=data_in';
elseif sz(2)==length(sym_list_in)
   
else
    error('incorrectly sized data')
end

data_aug=kron(data_in,ones(size(dd,1),1));
dd_aug=repmat(dd,[size(data_in,1),1]);
data=prod(data_aug.^dd_aug,2);
data=reshape(data,[size(dd,1),size(data_in,1)])';
end
