function [data]=monomials_data(data_in,degree)

%input list of data n_dataxn_terms, max total degree

%output: n_data row vectors off all monomial combinations of the input list from total degree 1 up to
%degree degree

%the function mss_asd is borrow from the Spotless toolbox
%https://github.com/spot-toolbox/spotless

nx=size(data_in,2);

dd=mss_asd(nx,1:degree);

data_aug=kron(data_in,ones(size(dd,1),1));
dd_aug=repmat(dd,[size(data_in,1),1]);
data=prod(data_aug.^dd_aug,2);
data=reshape(data,[size(dd,1),size(data_in,1)])';
end
