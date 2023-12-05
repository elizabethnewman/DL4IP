
clear; clc;


currDir = pwd;
dirName = '/Users/elizabethnewman/Library/CloudStorage/OneDrive-EmoryUniversity/[000] Code/[000] Public Libraries/IRtools';


cd(dirName)
addpath(genpath('.'))

cd(currDir);

currDir = pwd;
dirName = '/Users/elizabethnewman/Library/CloudStorage/OneDrive-EmoryUniversity/[000] Code/[000] Public Libraries/AIRtoolsII';


cd(dirName)
addpath(genpath('.'))

cd(currDir);

% options = {};
% count = 1;
% for tol = [1e-4,1e-6]
%     options{count} = DL4IP_OptionsDefault();
%     options{count}.VPAL.tol = tol;
%     count = count + 1;
% end


% % can be parallelized
% for i = 1:length(options)
%     [results,options,x,b] = DL4IP_RunExperiment(options);
% end

% setup problem
options           = DL4IP_OptionsDefault();
options.data.name = 'peppers';
options.data.size = [128,128];

options.dictionary.filename = '2023_10_17--flower_dictionary.mat';


options.AOperator.name   = 'deblur';
options.AOperator.params = {'downSamplingRatio', 4};
options.LOperator.name   = 'finite difference';
options.LOperator.lambda = 1e-2;

options.optimizer            = {'VPAL'};
options.VPAL.options.tol     = 1e-8;
options.VPAL.options.mu      = 1e-2;
options.VPAL.options.maxIter = 1000;
options.VPAL.options.display = 'iter';


% options.optimizer                   = {'MRNSDSparsity'};
% options.MRNSDSparsity.options.tol   = 1e-4;
% options.MRNSDSparsity.lambda        = 1e-8;
% options.MRNSDSparsity.options.MaxIter  = 200;

% options.optimizer            = {'MRNSDLight'};
% options.MRNSDLight.lambda     = 1e-8;
% options.MRNSDLight.options.tol     = 1e-8;
% options.MRNSDLight.options.mu      = 1e-3;
% options.MRNSDLight.options.maxIter = 1000;
% options.MRNSDLight.options.display = 'iter';


[results,options,x,b] = DL4IP_RunExperiment(options);

%% visualizations

XTrue = reshape(x,results{1}.probInfo.sizeX);
for i = 1:length(results)

    fig = figure(i); clf; fig.Name = results{i}.name;
    if strcmp(options.AOperator.name,'indicator')
        subplot(1,4,1);
        imagesc(results{1}.probInfo.B_dead);
    elseif strcmp(options.AOperator.name,'superresolution')
        figure(2);
        numberOfLowRes  = length(results{1}.probInfo.btrue);
        lowResImageSize = results{1}.probInfo.lowResImageSize;
        J = ceil(sqrt(numberOfLowRes));
        for j = 1:numberOfLowRes
            subplot(J,J,j)
            imshow(reshape(results{1}.probInfo.btrue{j},lowResImageSize(1),lowResImageSize(2)),[])
        end
         %  figure(1)
    else
        imagesc(reshape(b, results{1}.probInfo.sizeB));
    end
    
    % end
    ylabel('$\mathbf{b}$','Interpreter','latex')
   
    subplot(1,4,2)
    imagesc(XTrue);
    colorbar;
    cax = clim();
    ylabel('$\mathbf{x}$','Interpreter','latex')

    subplot(1,4,3);
    colorbar;
    imagesc(results{i}.XHat);
    clim(cax);
    colorbar;
    title(sprintf('rel err = %0.2e',results{i}.relErr))
    ylabel('$\widehat \mathbf{x}$','Interpreter','latex')

    subplot(1,4,4);
    imagesc(abs(results{i}.XHat - XTrue));
    colorbar;
    title(sprintf('nnz(X) / nnz(B) = %0.2f',results{i}.relNnzAlpha))
    ylabel('$\widehat \mathbf{x} - \mathbf{x}$','Interpreter','latex')
   
end
set(gca,'FontSize',18)
fig = figure(length(results) + 2); fig.Name = 'Rnrm';
for i = 1:length(results)
    semilogy(0:length(results{i}.optInfo.Rnrm),[results{i}.optInfo.initialEval.Rnrm;results{i}.optInfo.Rnrm(:)], 'LineWidth', 4, 'DisplayName', results{i}.name)
    hold on;
end
xlabel('iter')
ylabel('$||D * X - B||_F / ||B||_F$','Interpreter','latex')
legend()

hold off;

fig = figure(length(results) + 3); fig.Name = 'nnz';
for i = 1:length(results)
    plot(0:length(results{i}.optInfo.nnzX),[results{i}.optInfo.initialEval.nnzX;results{i}.optInfo.nnzX(:)] / nnz(XTrue), 'LineWidth', 4, 'DisplayName', results{i}.name)
    hold on;
end
xlabel('iter')
ylabel('nnz(X) / nnz(B)')
legend()
set(gca,'FontSize',18)
hold off;


fig = figure(length(results) + 4); fig.Name = 'alphas';
for i = 1:length(results)
    semilogy(1:length(results{i}.optInfo.alphas),results{i}.optInfo.alphas(:), 'LineWidth', 4, 'DisplayName', results{i}.name)
    hold on;
end
xlabel('iter')
ylabel('alpha')
legend()
set(gca,'FontSize',18)
hold off;



