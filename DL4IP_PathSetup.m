


if exist(dfile, 'file') ; delete(dfile); end

if ~exist('IRtools','dir')
    !git clone https://github.com/jnagy1/IRtools.git
end

if ~exist('AIRToolsII','dir')
    !git clone https://github.com/jakobsj/AIRToolsII
end


addpath(genpath('.'))
