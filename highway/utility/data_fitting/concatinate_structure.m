function data_out=concatinate_structure(data)

fnames=fieldnames(data);
for i=1:length(fnames)
    data_out.(fnames{i})=[];
end

for i=1:length(fnames)
    for j=1:size(data,1)
        for k=1:size(data,2)
           data_out.(fnames{i})=[data(j,k).(fnames{i});data(j,k).(fnames{i})];
        end
    end
end


end