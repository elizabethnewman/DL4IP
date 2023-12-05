
function[A,b,x,info] = deblurringProblem(X,varargin)

% band = bandwidth (size of blurring region)
% sigma = spread of blur (larger = more blurring)
% noiseLevel = amount of noise added after blurring

% set defaults
band        = 3;
sigma       = 4;  
noiseLevel  = 1e-4;   
plotting    = false;

% overwrite defaults
for k = 1:2:length(varargin)
    eval([varargin{k}, ' = varargin{k + 1};'])
end

% normalize
X = normalizeData(X,'Frobenius');
x = X(:);

% blurring operator
n  = size(X,1);
p  = [exp(-((0:band-1)) / (2 * sigma^2)),zeros(1,n-band)];
A1 = toeplitz(sparse(p));
A  = (1 / (2 * pi * sigma^2)) * kron(A1,A1);

b = A * X(:);
b = b + noiseLevel * norm(b(:)) * randn(size(b));

% problem info
info.sizeX = size(X);
info.sizeB = size(X);

if plotting
    % provide plots
    subplot(1,2,1)
    imshow(X,[]), title('true image')
    xlabel(num2str(n))
    ylabel(num2str(m))
    
   
    subplot(1,2,2)
    imshow(reshape(b,size(X)),[]), title('blurry, noisy image')
    drawnow
end

end
