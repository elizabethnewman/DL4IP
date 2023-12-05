function[A,b,x,info] = denoiseProblem(X,varargin)

% set defaults
noiseLevel = 1e-2;

% overwrite defaults
for k = 1:2:length(varargin)
    eval([varargin{k}, ' = varargin{k + 1};'])
end

% normalize and vectorize
X = normalizeData(X,'Frobenius');
x = X(:);

% setup operator
A = dOperator('identity',[size(X)]);
b = x + noiseLevel * norm(X(:)) * randn(size(X(:)));

info.sizeX = size(X);
info.sizeB = size(X);
end
