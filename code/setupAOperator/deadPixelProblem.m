function [A, b, x, info] = deadPixelProblem(X,varargin)

% set defaults
sparsity = 0.6;
plotting = false;

% overwrite
for k = 1:2:length(varargin)
    eval([varargin{k}, ' = varargin{k + 1};'])
end

% normalize
X = normalizeData(X,'max');
x = X(:);

% get image size
[m,n]   = size(X);      

l   = m*n;
idx = datasample(1:l,round(l*sparsity),'Replace',false); % get indices with dead pixel
B_dead      = X; 
B_dead(idx)  = 0;

% create operator with dead pixels
A        = speye(m*n); 
A(idx,:) = []; % create forward removal operator
b        = A * x;

info.sizeX      = size(X);
info.sizeB      = [size(A,1),size(X,2)];
info.deadPixels = idx; 
info.B_dead     = B_dead;

if plotting
    % provide plots
    subplot(1,2,1)
    imshow(X,[]), title(['true image (sparsity ',num2str(100*sparsity),'\%)'])
    xlabel(num2str(n))
    ylabel(num2str(m))
    
   
    subplot(1,2,2)
    imshow(B_dead,[]), title('observation (grey)')
    xlabel(['sparsity ',num2str(100*sparsity),'\%'])
    drawnow
end

end