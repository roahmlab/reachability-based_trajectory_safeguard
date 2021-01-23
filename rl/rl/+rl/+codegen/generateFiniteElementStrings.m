function [elementstr,IL,IR] = generateFiniteElementStrings(spec)
% GENERATEFINITEELEMENTSTRINGS
% generate an element string for finite set specs. Also return the
% appropriate brackets (IL,IR) if elements is a cell array or matrix

% Revised: 2-8-2019
% Copyright 2019 The MathWorks, Inc.


els = spec.Elements;
if iscell(els)
    % support elements defined by a cell array
    c = cellfun(@mat2str,els,'Uniformoutput',false);
    elementstr = sprintf('%s,',c{:});
    elementstr(end) = '';
    elementstr = ['{',elementstr,'}'];
    IL = '{';
    IR = '}';
else
    elementstr = mat2str(els);
    IL = '(';
    IR = ')';
end

