function[A,b,x,info] = superResolutionProblem(X_full,varargin)

% set defaults
downSamplingRatio   = 5;
numberOfLowRes      = 10;
plotting            = false;

% overwrite defaults
for k = 1:2:length(varargin)
    eval([varargin{k}, ' = varargin{k + 1};'])
end


[m,n] = size(X_full);                   % get size of full high image
X_full = X_full/max(X_full(:));           % normalize full high resolution image
xtrue = X_full(:);

% construct low resolution images
lowResImageSize = [m/downSamplingRatio, n/downSamplingRatio]; % number of blocks (including one reference image)

blockSize = lowResImageSize(1)*lowResImageSize(2);
mm = numberOfLowRes*blockSize;

HRparams.m1 = m;   HRparams.m2 = n;   % size of high resolution image
HRparams.h1 = 1/m; HRparams.h2 = 1/n; % pixel width
HRparams.ssp = downSamplingRatio;  % sub-sampling parameter

% creating non-linear parameter vector
t = linspace(-0.5, 0.5, numberOfLowRes); % rotation parameter
for i = 1:length(t)
  p1  = [0; t(i); 0];
  p2  = [-t(i);0;0];
  y{i} = [p1; p2];
end

% creating low resolution images and A block function handle
fprintf('Creating true low resolution images... ')
v = HRparams.h1; w = HRparams.h2;
[xx,yy] = ndgrid( v/2:v:(1-v/2) , w/2:w:(1-w/2) ); e = xx*0+1; PP = [xx(:),yy(:),e(:)];
fullbtrue = []; A = []; b = [];
for i = 1:length(t)
%   fprintf('%d, ',i)
  btrue{i} = getLRimage(xtrue,y{i},PP,HRparams);
  fullbtrue =  [fullbtrue; btrue{i}];
  Ablock = getAblock(y{i},HRparams);
  A = [A;Ablock];
  b = [b;btrue{i}];
end

fprintf(' done.\n')

if plotting
  fprintf('Creating plots...                      ')
  imshow(reshape(xtrue,m,n),[]), title('true superresolution image')
  J = ceil(sqrt(numberOfLowRes));
  figure
  for j = 1:numberOfLowRes
    subplot(J,J,j)
    imshow(reshape(btrue{j},lowResImageSize(1),lowResImageSize(2)),[])
  end
  title('low resolution image set')
  %   sgtitle('low resolution image set')
  drawnow
  fprintf(' done.\n')
end

x = xtrue;
info.X_full = X_full;
info.btrue  = btrue;
info.lowResImageSize = lowResImageSize;

end

function b = getLRimage(xtrue,y,PP,HRparams)
p1 = y(1:3); p2 = y(4:6);
u1 = PP*p1; u2 = PP*p2;
A = GetMRmat(u1,u2,HRparams);
b = A*xtrue;
end

function A = getAblock(y,HRparams)
v = HRparams.h1; w = HRparams.h2;
[xx,yy] = ndgrid( v/2:v:(1-v/2) , w/2:w:(1-w/2) ); e = xx*0+1; PP = [xx(:),yy(:),e(:)];
p1 = y(1:3); p2 = y(4:6);
u1 = PP*p1; u2 = PP*p2;
A = GetMRmat(u1,u2,HRparams);
end


function[A,S,R] = GetMRmat(u1,u2,param)
% [A,S,R] = GetMRmat(u1,u2,param)
%
% This function gets the super-resolution matrix corresponding to given
% displacements.
%
%   Input: u1, u2 - displacements in the x and y directions
%           param - parameters
%
%   Output:     A - MR matrix (R*S)
%               S - bilinear interpolation matrix of weights
%               R - restriction matrix
%

h1 = param.h1; h2 = param.h2; m1 = param.m1; m2 = param.m2;
xnodal = 0:h1:m1*h1; xmid =  xnodal(1:end-1) + diff(xnodal)/2;
ynodal = 0:h2:m2*h2; ymid =  ynodal(1:end-1) + diff(ynodal)/2;


[X,Y] = ndgrid(xmid,ymid);

S = Bilinear2(m1,m2,X+reshape(u1,m1,m2),Y+reshape(u2,m1,m2),h1,h2);

R = rest2D(m1,m2,param.ssp);
A = R*S;

end

function A = Bilinear2(m,n,x,y,h1,h2)
% A = Bilinear2(m,n,x,y,h1,h2)
%
% This functions uses Bilinear interpolation to produce a sparse matrix of
% "weights"
%
% ***********************************************************************
%   The difference between Bilinear2 and Bilinear is that Bilinear2 takes
%   into account the boundary and does NOT assume that the image is
%   embedded in a zero boundary.
% ***********************************************************************
%   Input:
%        m,n - size of the HR image
%        x,y - displaced pixels
%      h1,h2 - pixel size
%
%   Output:
%          A - sparse matrix of coefficients

% Convert x and y to the coordinate system 1:m, 1:n
x = x/h1 + 1/2; y = y/h2 + 1/2;

% Vectorized version
%     Th = zeros(m*n,1);
j=floor(x(:)); xi  = x(:)-j;
k=floor(y(:)); eta = y(:)-k;

ind1 = find(1<=j & j<m & 1<=k & k<n);
jk = j(ind1) + (k(ind1)-1)*m;
j1k = j(ind1) + 1+ (k(ind1)-1)*m;
jk1 = j(ind1) + k(ind1)*m;
j1k1 = j(ind1) + 1 + k(ind1)*m;

ii = [ind1;ind1;ind1;ind1];
jj = [jk;j1k;jk1;j1k1];
ss = [(1-xi(ind1)).*(1-eta(ind1)); (xi(ind1)).*(1-eta(ind1)); ...
  (1-xi(ind1)).*(eta(ind1));   (xi(ind1)).*(eta(ind1))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind2 = find(j==m & 1<k & k<n);
jk = j(ind2) + (k(ind2)-1)*m;
jk1 = j(ind2) + k(ind2)*m;

ii = [ii; ind2;ind2];
jj = [jj; jk;jk1];
ss = [ss; (1-eta(ind2)); (eta(ind2))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind3 = find(1<j & j<m & k==n);
jk = j(ind3) + (k(ind3)-1)*m;
j1k = j(ind3) + 1+ (k(ind3)-1)*m;

ii = [ii; ind3;ind3];
jj = [jj; jk;j1k];
ss = [ss; (1-xi(ind3)); (xi(ind3))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind4 = find(j==m & k==n);
jk = j(ind4) + (k(ind4)-1)*m;

ii = [ii; ind4];
jj = [jj; jk];
ss = [ss; ind4./ind4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind5 = find(j==m & k==1);
jk = j(ind5) + (k(ind5)-1)*m;

ii = [ii; ind5];
jj = [jj; jk];
ss = [ss; ind5./ind5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind6 = find(j==1 & k==n);
jk = j(ind6) + (k(ind6)-1)*m;

ii = [ii; ind6];
jj = [jj; jk];

ss = [ss; ind6./ind6];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = sparse(ii,jj,ss,n*m,n*m);

end

function [R] = rest2D(m,n,ssp)
%
%   [R] = rest2D(m,n,ssp)
%
%   This function computes the restriction operator for the
%       Super-Resolution codes.
%
%

if ~exist('ssp','var'), ssp = 2; end

if rem(n,ssp)~=0
  error('N must divide sub sampling');
end

D1 = kron(speye(n/ssp),ones(1,ssp));
D2 = kron(speye(m/ssp),ones(1,ssp));

R = kron(D1,D2)/ssp^2;
end



