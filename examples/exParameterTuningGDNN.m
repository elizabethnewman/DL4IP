
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

options.GDNN.options.verbose = false;
options.GDNN.options.MaxIter = 100;

RESULTS = {};
count   = 0;
aRange  = 0.1:0.1:1;
cRange  = -1:0.1:0;
for a = aRange
    for c = cRange
        fprintf('count = %d: Starting a = %0.4f, c = %0.4f...',count,a,c)
        w  = @(x) (x - c > 0) .* (exp(a * (x - c)) - a * (x - c) - 1);
        dw = @(x) (x - c > 0) .* (a * exp(a * (x - c)) -  a);

        startTime       = tic; 
        [coeff,optInfo] = gradientDescentPointwiseMapping(AD, [b(:);zeros(size(AD.L,1),1)],w,dw,options.GDNN.options);
        endTime         = toc(startTime);
        
        % compute results
        results = computeResults(coeff,optInfo,probInfo,x,AD.D,endTime,'GDNN',options);

        RESULTS{count + 1} = results;
        count              = count + 1;
        fprintf('...finished!\n')
    end
end

%% 
save([saveDirName,'/RESULTS'],'RESULTS','aRange','cRange')

return;

%% 

count  = 0;
RELERR = [];
RELNNZ = [];
for i = 1:length(aRange)
    for j = 1:length(cRange)
        RELERR(i,j) = RESULTS{count + 1}.relErr;
        RELNNZ(i,j) = RESULTS{count + 1}.relNnzAlpha;
        count = count + 1;
    end
end

%%

[A,C] = meshgrid(aRange,cRange);
[AA,CC] = meshgrid(...
    linspace(aRange(1)-0.05,aRange(end)+0.05,length(aRange)+1),...
    linspace(cRange(1)-0.05,cRange(end)+0.05,length(cRange)+1));

fig = figure(1); clf;
hi = imagesc(cRange(:),aRange(:),RELERR);
hold on;
hm = mesh(CC,AA,zeros(size(AA)));
hm.FaceColor = 'none';
hm.EdgeColor = 'k';
% colorbar;
clim([0,1])
% xlabel('c')
% ylabel('a')
axis('off')
axis('image')
set(gca,'YDir','normal')
hold off;
exportgraphics(fig,[saveDirName,'relerr.png'],'BackgroundColor','none')


figure(1); clf;
hi = imagesc(cRange(:),aRange(:),RELNNZ);
hold on;
hm = mesh(CC,AA,zeros(size(AA)));
hm.FaceColor = 'none';
hm.EdgeColor = 'k';
% colorbar;
clim([0,1])
% xlabel('c')
% ylabel('a')
axis('off')
axis('image')
set(gca,'YDir','normal')
hold off;
exportgraphics(fig,[saveDirName,'sparsity.png'],'BackgroundColor','none')



