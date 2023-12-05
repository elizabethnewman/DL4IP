function[b] = dictionaryOperator(alpha,transFlag,D,patchSize,imgSize)
%
% Let X be an M x N image
% Let Y be the patchified version of X size pq x (M/p)(N/q) (non-overlapping patches)
%
% We can write Y \approx D * Alpha where D is the dictionary and Alpha are
% the coefficients.
%
% In vectorized form, we have y \approx kron(I,D) * alpha, which is how we
% want to think about this function.
% 
% Inputs:
%   alpha       : vectorized coefficients (non-transposed), vectorized image (transposed)
%   transFlag   : boolean or string 'transp' indicating transpose operation or not (default=0)
%   D           : dictionary of size pq x s (no Kronecker product)
%   patchSize   : 1 x 2 array of patch size [p,q]
%   imgSize     : 1 x 2 array of image size [M,N]
%   
%

if nargin == 0, runMinimalExample; return; end

if exist('transFlag','var') && (strcmp(transFlag, 'transp') || all(transFlag == 1))
    b = reshape(D' * im2col(reshape(alpha,imgSize),patchSize,'distinct'),[],1);
else
    b = reshape(col2im(D * reshape(alpha,size(D,2),[]),patchSize,imgSize,'distinct'),[],1);
end

end


function runMinimalExample()

[p,q] = deal(10,10);
[M,N] = deal(50,50);

s     = 100;
D     = rand(p * q, s);
alpha = rand(s, floor(M / p) * floor(N / q));

b = dictionaryOperator(alpha(:), 0, D, [p,q], [M,N]);

alpha2 = dictionaryOperator(b(:), 0, D, [p,q], [M,N]);

end


