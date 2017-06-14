function loaddeps()
    [current, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(current, 'condense'));
    addpath(fullfile(current, 'cytof'));
    addpath(fullfile(current, 'phate/Matlab'));
    addpath(fullfile(current, 'sanky'));
    addpath(fullfile(current, 'sanky/template'));
end