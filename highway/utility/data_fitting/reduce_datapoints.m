function out_data=reduce_datapoints(in_data,reduce_num)
names=fieldnames(in_data);
out_data=struct();

for i=1:length(names)
    temp_data=getfield(in_data,names{i});
    V=temp_data(1:reduce_num:length(in_data.time));
    out_data=setfield(out_data,names{i},V);
end

end