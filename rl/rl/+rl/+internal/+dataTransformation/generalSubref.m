function Data = generalSubref(Data, SubrefIdx, DimToSlice)
%GENERALSUBREF indexes Data on dimension DimToSlice (usually Batch 
%dimension) with Idx.
%   Data is a cell array of numerical matrices.
%   DimToSlice is a cell array of scalar. DimToSlice{i} <= ndims(Data{i})
%   SubrefIdx is a numeric vector indicates the subref indexes.
%   Number of cell for Data and DimToSlice must match
%       e.g. if size(InputData{1}) = [x1 x2 x3 x4] => 
%       InputData = generalSubref(Data, SubrefIdx, {3}) is equivalent to
%           OutputData{1} = InputData{1}(:,:,SubrefIdx,:)
%           OutputData{2} = InputData{2}(:,:,SubrefIdx,:)
%           OutputData{n} = InputData{n}(:,:,SubrefIdx,:)

% Copyright 2019 The MathWorks Inc.

if iscell(Data)
    Data = cellfun(@(x,y) iSliceDataFromSubrefIdx(x, SubrefIdx, y),Data,DimToSlice,'UniformOutput',false);
else
    Data = iSliceDataFromSubrefIdx(Data, SubrefIdx, DimToSlice);
end

end

function Data = iSliceDataFromSubrefIdx(Data, Idx, DimToSlice)
LastDims = ndims(Data);
[SubsciptVal{1:LastDims}] = deal(':');
SubsciptVal{DimToSlice} = Idx;
Data = subsref(Data,substruct('()',SubsciptVal));
end
