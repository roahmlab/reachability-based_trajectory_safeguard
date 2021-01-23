function rlCreateEnvTemplate(className)
% rlCREATEENVTEMPLATE: Makes a copy of RL environment template class

% Copyright 2018 The MathWorks, Inc.

% Error checking
validateattributes(className,{'char','string'},{'nonempty','scalartext'});
if className == ""
    error(message('rl:general:expectedScalartext'));
end
className = char(className);
% Path to the MATLAB template class
filepath = fullfile('+rl','+env','TemplateForMATLABEnv.m');
% Read template and rename class and constructor
textBuffer = readMatfile(filepath);
textBuffer = regexprep(textBuffer,'TemplateForMATLABEnv',className);
textBuffer = regexprep(textBuffer,upper('TemplateForMATLABEnv'),upper(className));
% Open new custom class on editor
matlab.desktop.editor.newDocument(textBuffer);

function textBuffer = readMatfile(filepath)
    fid = fopen(filepath);
    textBuffer = '';
    while(~feof(fid))
        textBuffer = [textBuffer fgetl(fid) newline];
    end
    fclose(fid);
end

end