function data_out=filter_logic(data,L)
data_out=data;
fnames=fieldnames(data);
for i=1:length(fnames)
    valin=data.(fnames{i});
    if length(data.(fnames{i}))==length(L)
    data_out.(fnames{i})=valin(L,:);
    end
end


end