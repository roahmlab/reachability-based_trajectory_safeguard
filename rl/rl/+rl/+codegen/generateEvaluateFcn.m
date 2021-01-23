function varargout = generateEvaluateFcn(argStruct)
% GENERATEEVALUATEFCN
% generate an N function for creating a static function equivalent to
% calling predict on a representation

% Revised: 2-4-2019
% Copyright 2019 The MathWorks, Inc.

% parse the struct
fcnname         = argStruct.FunctionName;
inputstr        = argStruct.InputString;
outputstr       = argStruct.OutputString;
bodystr         = argStruct.BodyString;
localfcnstr     = argStruct.LocalFunctionString;
outputstr_y     = rl.codegen.handleMultipleOutputStrings(outputstr);

topstr = sprintf([...
    'function %s = %s(%s)\n',...
    '%%#codegen',...
    '\n\n',...
    '%% Reinforcement Learning Toolbox\n',...
    '%% Generated on: %s'],...
    outputstr_y,fcnname,inputstr,datestr(now));
middlestr = sprintf('%s\nend',bodystr);
bottomstr = sprintf('%s',localModLocalFcn(localfcnstr));
str = sprintf('%s\n\n%s\n%%%% Local Functions\n%s',topstr,middlestr,bottomstr);

% remove trailing carriage returns
while strcmp(str(end),newline)
    str(end) = '';
end

if nargout
    varargout{1} = str;
else
    % print to file
    filename = fcnname;
    if ~strcmp(filename(end:-1:end),'.m')
        filename = [filename,'.m'];
    end
    fid = fopen(filename,'w');
    fwrite(fid,str);
    fclose(fid);
end

function str = localModLocalFcn(str)
if ~isempty(str)
    str = sprintf('%s\n\n',str);
end
