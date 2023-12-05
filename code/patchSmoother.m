function[L] = patchSmoother(patchSize,imgSize)
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

if nargin == 0; runMinimalExample(); return; end

% get sizes
M  = imgSize(1);     N  = imgSize(2);
p  = patchSize(1);   q  = patchSize(2);
Mp = floor(M / p);   Nq = floor(N / q);

% number of patches
numPatch = Mp * Nq;


% smooth patches top-to-bottom (applied from right to patchified image)
D1               = sparse(Mp,Mp - 1);
idxDiagD         = (1:Mp:numel(D1)) + ((1:size(D1,2))-1);
D1(idxDiagD)     = 1;
D1(idxDiagD + 1) = -1;
L1               = kron(speye(Nq), D1);

% p * (q - 1) = number of patch boundaries across left-to-right
L2               = sparse(numPatch,Mp * (Nq - 1));
idxDiag          = (1:numPatch:numel(L2)) + ((1:size(L2,2))-1);
L2(idxDiag)      = 1;
L2(idxDiag + Mp) = -1;

L = [L1,L2];
% 
% 
% idx_q = q:q:N-q;
% % idx_q = sort([idx_q,idx_q+1],'ascend');
% % idx_p = kron([p,p+1],1:floor(M/p));
% idx_p = p:p:M-p;
% % idx_p = sort([idx_p,idx_p+1],'ascend');
% 
% 
% F1 = sparse(length(idx_p),M);
% F1(sub2ind([length(idx_p),M],1:length(idx_p),idx_p))     = 1;
% F1(sub2ind([length(idx_p),M],1:length(idx_p),idx_p + 1)) = -1;
% % F1 = F1(1:end-1,:);
% % F1 = F1(idx_p,:);
% 
% 
% % F2 = eye(N) - diag(ones(N-1,1),-1,N,N);
% % F2 = spdiags([ones(N,1),-ones(N,1)],[-1,0],N,N);
% % F2 = F2(:,1:end-1);
% % F2 = F2(:,idx_q);
% F2 = sparse(N,length(idx_q));
% F2(sub2ind([N,length(idx_q)],idx_q,1:length(idx_q)))     = 1;
% F2(sub2ind([N,length(idx_q)],idx_q+1,1:length(idx_q))) = -1;
% 
% 
% L = [kron(eye(N),F1); kron(F2',eye(M))];


end


function runMinimalExample()

% construct patch smoothing operator
p = 3;
q = 4;

M = 18;
N = 20;

L = patchSmoother([p,q],[M,N]);


n1 = floor(M / p);
n2 = floor(N / q);
A  = zeros(p * q, n1 * n2);

for i = 1:size(A,2)
    if mod(i,4)
        A(:,i) = i * kron((1:p)',ones(q,1));
    else
        A(:,i) = i * kron(ones(p,1),(1:q)');
    end
end
tmp = col2im(A,[p,q],[M,N],'distinct');
tmp2 = A * L;

F1 = eye(M) + diag(ones(M-1,1),1);
F1 = F1(1:end-1,:);

F2 = eye(N) + diag(ones(N-1,1),-1);
F2 = F2(:,1:end-1);

L_FD = [kron(eye(N),F1); kron(F2',eye(M))];

% patch
tt      = reshape(1:numel(A),M,[]);
tt2     = im2col(tt,[p,q],'distinct');
[~,idx] = sort(tt2(:),'ascend');
% gives norm(tmp(:) - A(idx)) = 0

% tt      = reshape(1:numel(A),size(A));
% tt2     = col2im(tt,[p,q],[M,N],'distinct');
% [~,idx] = sort(tt2(:),'ascend');
% $  gives norm(tmp(idx) - A(:)) = 0

L_P = kron(L',eye(size(A,1)));
L_P = L_P(:,idx);

norm(tmp2(:) - L_P * tmp(:))



figure(1); clf;
subplot(2,2,1);
imagesc(tmp); colorbar;
cax = clim();
title('A')

subplot(2,2,2);
imagesc(A); colorbar;
title('patchify(A)')

subplot(2,2,3);
spy(L);
title('L')

subplot(2,2,4);
imagesc(tmp2);
title('patchify(A) * L')
colorbar;



end


