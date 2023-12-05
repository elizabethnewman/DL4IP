function[results] = computeResults(alpha,optInfo,probInfo,x,D,endTime,name,options)

patchSize = probInfo.patchSize;
imgSize   = probInfo.sizeX;
X         = reshape(x,imgSize);

results.name            = name;
results.optInfo         = optInfo;
results.probInfo        = probInfo;
results.options         = options;
results.alpha           = reshape(alpha,size(D,2),[]); % coefficients
results.XHat            = col2im(D * results.alpha,patchSize,imgSize,'distinct');
results.err             = norm(X - results.XHat,'fro');
results.relErr          = results.err / norm(X,'fro');
results.nnzAlpha        = nnz(results.alpha);
results.relNnzAlpha     = results.nnzAlpha  / nnz(X);
results.time            = endTime;

end

