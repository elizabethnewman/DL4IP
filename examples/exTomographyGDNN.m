
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
options.data.name = 'tomography';
options.data.size = [256,256];

% choose dictionary
options.dictionary.filename = '2023_10_17--flower_dictionary.mat';
% options.dictionary.filename = '08-29-23-ScheppLoganDictionary_ADMM_patch.mat';

% operator
options.AOperator.name   = 'tomography';
options.LOperator.name   = 'patch smoother';
options.LOperator.lambda = 1e0;

% set seed for reproducibility
rng(options.seed)

% setup problem
[AD,x,b,probInfo] = DL4IP_ProblemSetup(options);

%% nonnegative mapping
[a,c]   = deal(0.5,-1);
w       = @(x) (x - c > 0) .* (exp(a * (x - c)) - a * (x - c) - 1);
dw      = @(x) (x - c > 0) .* (a * exp(a * (x - c)) -  a);

options.GDNN.options.verbose = true;
options.GDNN.options.MaxIter = 100;

startTime       = tic; 
[coeff,optInfo] = gradientDescentPointwiseMapping(AD, [b(:);zeros(size(AD.L,1),1)],w,dw,options.GDNN.options);
endTime         = toc(startTime);

% compute results
results = computeResults(coeff,optInfo,probInfo,x,AD.D,endTime,'GDNN',options);

%% save results

save([saveDirName,'/results'],'results')


% stopping criteria
numSC  = size(optInfo.values,2) - length(optInfo.headers) + 1;
optInfo.headers{end} = 'SC1';
for i = 2:numSC
    optInfo.headers = cat(2,optInfo.headers,{['SC',num2str(i)]});
end

T = array2table(optInfo.values,'VariableNames',cat(2,optInfo.headers));
writetable(T,[saveDirName,'results.csv']);

%% plot approximations

set(0,'DefaultFigureWindowStyle','normal')

fig = figure(1); clf;
imagesc(reshape(x,results.probInfo.imgSize)); colormap parula;
cax = clim();
axis('off')
axis('image')
exportgraphics(fig,[saveDirName,'orig.png'],'BackgroundColor','none')

fig = figure(1); clf;
imagesc(reshape(results.XHat,results.probInfo.imgSize)); colormap parula;
clim(cax);
axis('off')
axis('image')
exportgraphics(fig,[saveDirName,'recon.png'],'BackgroundColor','none')

fig = figure(1); clf;
imagesc(reshape(abs(results.XHat(:) - x),results.probInfo.imgSize)); colormap hot;
% clim([0,1e1]);
clim([0,1.458088334803294])
% colorbar;
axis('off')
axis('image')
exportgraphics(fig,[saveDirName,'abs_diff.png'],'BackgroundColor','none')

%% plot convergence

fig = figure(1); clf;
semilogy(optInfo.values(:,1),optInfo.values(:,4),'-o','LineWidth',3)
xlabel('iter')
ylabel('f')
set(gca,'FontSize',18)
matlab2tikz([saveDirName,'Rnrm.tex'])


fig = figure(1); clf;
semilogy(optInfo.values(:,1),optInfo.values(:,10),'-o','LineWidth',3)
xlabel('iter')
ylabel('alpha')
set(gca,'FontSize',18)
matlab2tikz([saveDirName,'alphas.tex'])

%% store solution and operator
fig = figure(1); clf;
spy(results.alpha)
axis('off')
exportgraphics(fig,[saveDirName,'solution.png'],'BackgroundColor','none')


fig = figure(1); clf;
spy(AD.A)
axis('off')
exportgraphics(fig,[saveDirName,'A.png'],'BackgroundColor','none')


fig = figure(1); clf;
spy(AD.L)
axis('off')
exportgraphics(fig,[saveDirName,'A.png'],'BackgroundColor','none')


%% 
diary off;

