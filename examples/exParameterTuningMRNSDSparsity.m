
clear; clc;


saveDirName = [mfilename,'--',datestr(datetime('today')),'/'];
if ~exist(saveDirName,'dir'), mkdir(saveDirName); end


copyfile([mfilename,'.m'],saveDirName);

dfile = [saveDirName,mfilename,'.txt'];
if exist(dfile, 'file') ; delete(dfile); end
diary(dfile)

%% setup paths

run('../DL4IP_PathSetup')

%%
% setup problem
options           = DL4IP_OptionsDefault();
options.data.name = 'peppers';
options.data.size = [256,256];

% choose dictionary
options.dictionary.filename = '2023_10_17--flower_dictionary.mat';

% operator
options.AOperator.name   = 'deblur';
options.LOperator.name   = 'finite difference';
options.LOperator.lambda = 1e-4;


% set seed for reproducibility
rng(options.seed)

% setup problem
[AD,x,b,probInfo] = DL4IP_ProblemSetup(options);

%% nonnegative mapping
[a,c]   = deal(0.1,-0.75);
w       = @(x) (x - c > 0) .* (exp(a * (x - c)) - a * (x - c) - 1);
dw      = @(x) (x - c > 0) .* (a * exp(a * (x - c)) -  a);

options.MRNSDSparsity.options.MaxIter = 100;
options.MRNSDSparsity.options.x0      = w(rand(size(AD,2),1));



RESULTS = {};
count   = 0;
lambdaRange  = sort([0,logspace(-12,-4,17)]);
for lambda = lambdaRange
    fprintf('count = %d: Starting lambda = %0.4e...',count,lambda)
    options.MRNSDSparsity.lambda = lambda;

    startTime       = tic; 
    [coeff,optInfo] = IRmrnsdSoftThreshold(AD, [b(:);zeros(size(AD.L,1),1)],options.MRNSDSparsity.lambda, options.MRNSDSparsity.options);
    endTime         = toc(startTime);
    
    % compute results
    results = computeResults(coeff,optInfo,probInfo,x,AD.D,endTime,'GDNN',options);

    RESULTS{count + 1} = results;
    count              = count + 1;
    fprintf('...finished!\n')
end

%% 
save([saveDirName,'/RESULTS'],'RESULTS','lambdaRange')

return;

%% 

count  = 0;
RELERR = [];
RELNNZ = [];
for i = 1:length(lambdaRange)
    RELERR(i) = RESULTS{count + 1}.relErr;
    RELNNZ(i) = RESULTS{count + 1}.relNnzAlpha;
    count = count + 1;
end

T = array2table([lambdaRange(:),RELERR(:),RELNNZ(:)],'VariableNames',{'lambda','rel_err','rel_sparsity'});
writetable(T,[saveDirName,'RESULTS.csv']);


%%
fig = figure(1); clf;

semilogx(lambdaRange,RELERR)
hold on;
semilogx(lambdaRange,RELNNZ)
% yline(0.1771,'b--')
% yline(0.4009,'r--')

hold off;



% semilogx(lambdaRange(8:14),RELERR(8:14))
% hold on;
% semilogx(lambdaRange(8:14),RELNNZ(8:14))
% hold off;
% matlab2tikz([saveDirName,'Rnrm.tex'])

% exportgraphics(fig,[saveDirName,'relerr.png'],'BackgroundColor','none')




