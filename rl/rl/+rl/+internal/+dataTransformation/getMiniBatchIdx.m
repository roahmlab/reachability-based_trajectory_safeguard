function MiniBatchIdx = getMiniBatchIdx(FullLength, MiniBatchSize, IsShuffle)
%GETMINIBATCHIDX returns a cell array of indices. Each cell contains 
%   MiniBatchSize indices. The total number of indices in all cell is FullLength
%   IsShuffle determines if the indices are shuffled, else returns sorted
%   indices. Assume MiniBatchSize > 0. Return empty matrix if FullLength <= 0
%   e.g. MiniBatchIdx = rl.util.getMiniBatchIdx(5,2,false)

% Copyright 2019 The MathWorks Inc.

if nargin < 3
    IsShuffle = true;
end

if IsShuffle
    Idx = randperm(FullLength);
else
    Idx = 1:FullLength;
end

NumBatch = ceil(FullLength / MiniBatchSize);
MiniBatchIdx = cell(NumBatch, 1);
for ct = 1:NumBatch
    StartIdx = (ct - 1) * MiniBatchSize + 1;
    EndIdx   = min(ct * MiniBatchSize, FullLength);
    MiniBatchIdx{ct} = Idx(StartIdx:EndIdx);
end

end

