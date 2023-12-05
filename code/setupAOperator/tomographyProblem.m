function[A,b,x,info] = tomographyProblem(X,varargin)
% requires AIRToolsII

% set defaults
n_projections  = 100;
noiseLevel     = 1e-4;            

% overwrite defaults
for k = 1:2:length(varargin)
    eval([varargin{k}, ' = varargin{k + 1};'])
end

% create data (no normalization)
x = X(:);

% create operator
angles           = linspace(0,179,n_projections);
[A,b,~,~,n_rays] = paralleltomo(size(X,1),angles); % load tomo operator

% add noise
b = b + noiseLevel * norm(b(:)) * randn(size(b));

% store info
info.sizeX = size(X);
info.sizeB = [n_rays,n_projections];
end