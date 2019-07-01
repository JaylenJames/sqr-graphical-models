function [ baseFilename, baseDir ] = createfilename(args, argValues, mfilenameStr)
%CREATEFILENAME Determine base filename and directory from executing function file
%and arguments to the function (for running experiments)

% Determine base directory
[funcDir, funcFilename] = fileparts(mfilenameStr);
if(ispc() || ismac()) 
    gitInfo = 'GIT_UNKNOWN';
else
    [~, gitDate] = system('export TERM=ansi; git log -n 1 --date=iso --format="%cd"');
    gitDate = strtrim(gitDate);
    gitDate = strrep(strrep(gitDate(1:16),':','-'), ' ', '_');
    [~, gitHash] = system('export TERM=ansi; git log -n 1 --date=iso --format="%h"');
    gitHash = strtrim(gitHash);
    gitHash = gitHash(1:7);
    gitInfo = sprintf('%s_%s', gitDate, gitHash); % date and hash
end
%topDir = fileparts(funcDir); % Assume executed in code directory
baseDir = fullfile('..', 'data', 'generated', sprintf('%s_%s', funcFilename, gitInfo));

% Determine filename
filename = funcFilename;
for i = 1:size(args,1)
    argStruct = args(i,1);
    value = argValues{i};
    filename = sprintf('%s_%s', filename, arg2string(argStruct, value));
end
baseFilename = sprintf('%s_%s', filename, gitInfo);

end

function s = arg2string(argStruct, value)
    valString = 'UNDEFINED';
    if(isnumeric(value))
        valString = sprintf('%g', value);
    elseif(ischar(value))
        valString = value;
    elseif(islogical(value))
        if(value); valString = 'T';
        else valString = 'F'; end;
    end
    NUM_CHAR = 4;
    if(length(argStruct.name) <= NUM_CHAR)
        nameString = argStruct.name;
    else
        nameString = argStruct.name(1:NUM_CHAR);
    end
    s = sprintf('%s-%s', nameString, valString);
end


