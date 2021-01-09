function data=filter_field(data,datafield,range)
L=datafield>=range(1)&datafield<=range(2);
fnames=fieldnames(data);
for i=1:length(fnames)
    valin=data.(fnames{i});
    data.(fnames{i})=valin(L,:);
end


end