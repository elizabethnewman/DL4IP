
function[options] = DL4IP_OptionsDefault()

% reproducibility
options.seed       = 42;

% dictionary
% 2023_10_17--flower_dictionary.mat
% 2023_10_24--earth_dictionary
options.dictionary.filename = '2023_10_17--flower_dictionary.mat';

% data
% {'moon', 'GirlWithPearlEarring', 'peppers', 'tomography'}
options.data.name       = 'moon';
options.data.size       = [128,128];
options.data.color      = false;
options.data.seed       = 123;

% forward operator
% {'deblur', 'indicator', 'superresolution', 'denoise', 'tomography'}
options.AOperator.name   = 'deblur';
options.AOperator.params = {};

% regularization operator 
% {'none', 'Tikhonov', 'finite difference', 'patch smoother'}
options.LOperator.name    = 'none';
options.LOperator.lambda  = 1e-2;

% optimizers
options.optimizer = 'MRNSDL1';

options.MRNSDL1.options = struct('x0','none', ...
    'MaxIter',100, 'x_true','none', ...
    'NoiseLevel','none', 'eta',1.01, ...
    'NE_Rtol',1e-12, 'IterBar','on', ...
    'NoStop', 'off');
options.MRNSDL1.lambda  = 2e-5;  % sparsity regularization


options.MRNSDSparsity.options = struct('x0','none', ...
    'MaxIter',100, 'x_true','none', ...
    'NoiseLevel','none', 'eta',1.01, ...
    'NE_Rtol',1e-12, 'IterBar','on', ...
    'NoStop', 'off');
options.MRNSDSparsity.lambda = 1e-8;  % sparsity regularization


options.VPAL.options = struct('D',[],'maxIter',[], 'display', 'iter',...
    'tol', 1e-6, 'xtrue', [], 'lambda', 1, 'mu', 0, 'getAllInfo',0,...
    'stepSize','linearized','dof',0,'bnd',[0,inf]);

options.GDNN.options = struct('z0','none', 'MaxIter',100, 'x_true','none', ...
    'NoiseLevel','none', 'eta',1.01, 'NE_Rtol',1e-12, ...
    'absTol', 1e-8, 'relTol', 1e-8,...
    'NoStop', 'off','verbose',false);

end



