function [data,b]=filterout_nans_and_infs(data,b)
dataT=data';
nan_rows=any(isnan(dataT));
data=data(~nan_rows,:);
b=b(~nan_rows);
dataT=data';
inf_rows=any(abs(dataT)==Inf);
data=data(~inf_rows,:);
b=b(~inf_rows);
end