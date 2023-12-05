

function[results,options,x,b] = DL4IP_RunExperiment(options)

% set seed for reproducibility
rng(options.seed)

% setup problem
[AD,x,b,probInfo] = DL4IP_ProblemSetup(options);

count = 1;

% run optimizers
if any(contains(options.optimizer,'MRNSDL1'))
    % options.MRNSDL1.options.MaxIter = 10 * size(b,1);

    startTime       = tic; 
    [alpha,optInfo] = IRmrnsdL1(AD, [b(:);zeros(size(AD.L,1),1)], options.MRNSDL1.lambda, options.MRNSDL1.options);
    endTime         = toc(startTime);
    results{count}  = computeResults(alpha,optInfo,probInfo,x,AD.D,endTime,'MRNSDL1');
    count = count + 1;
end

if any(contains(options.optimizer,'MRNSDSparsity'))    
    % options.MRNSDSparsity.options.MaxIter = 10 * size(b,1);
    
    startTime       = tic; 
    [alpha,optInfo] = IRmrnsdSoftThreshold(AD, [b(:);zeros(size(AD.L,1),1)], options.MRNSDSparsity.lambda, options.MRNSDSparsity.options);
    endTime         = toc(startTime);
    results{count}  = computeResults(alpha,optInfo,probInfo,x,AD.D,endTime,'MRNSDSparsity');
    count = count + 1;
end


if any(contains(options.optimizer,'MRNSDLight'))    
    % options.MRNSDSparsity.options.MaxIter = 10 * size(b,1);
    
    startTime       = tic; 
    [alpha,optInfo] = mrnsdLight(AD, [b(:);zeros(size(AD.L,1),1)], options.MRNSDLight.lambda, options.MRNSDLight.options);
    endTime         = toc(startTime);
    results{count}  = computeResults(alpha,optInfo,probInfo,x,AD.D,endTime,'MRNSDLight');
    count = count + 1;
end

if any(contains(options.optimizer,'VPAL'))
    D                            = dOperator('identity',size(AD,2));
    options.VPAL.options.D       = D;
    
    if isempty(options.VPAL.options.maxIter)
        options.VPAL.options.maxIter = 10 * size(b,1);
    end
    
    startTime         = tic; 
    [alpha,~,optInfo] = vpalmrnsd(AD, [b(:);zeros(size(AD.L,1),1)],struct2cellWithFields(options.VPAL.options));
    endTime           = toc(startTime);
    results{count}    = computeResults(alpha,optInfo,probInfo,x,AD.D,endTime,'VPAL');
end


end

%% helper functions
function[results] = computeResults(alpha,optInfo,probInfo,x,D,endTime,name)

patchSize = probInfo.patchSize;
imgSize   = probInfo.sizeX;
X         = reshape(x,imgSize);

results.name            = name;
results.optInfo         = optInfo;
results.probInfo        = probInfo;
% results.options         = options;
results.alpha           = reshape(alpha,size(D,2),[]); % coefficients
results.XHat            = col2im(D * results.alpha,patchSize,imgSize,'distinct');
results.err             = norm(X - results.XHat,'fro');
results.relErr          = results.err / norm(X,'fro');
results.nnzAlpha        = nnz(results.alpha);
results.relNnzAlpha     = results.nnzAlpha  / nnz(X);
results.time            = endTime;

end

function[C] = struct2cellWithFields(S)

C = {};

F = fields(S);
for i = 1:length(F)
    C = cat(2,C,F{i});
    C = cat(2,C,{S.(F{i})});
end
end

