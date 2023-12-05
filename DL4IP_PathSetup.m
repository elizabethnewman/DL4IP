


if exist(dfile, 'file') ; delete(dfile); end

if ~exist('IRtools','dir')
    !git clone https://github.com/jnagy1/IRtools.git
end

if ~exist('AIRToolsII','dir')
    !git clone https://github.com/jakobsj/AIRToolsII
end

if ~exist('matlab2tikz','dir')
    !git clone https://github.com/matlab2tikz/matlab2tikz
end


addpath(genpath('.'))
