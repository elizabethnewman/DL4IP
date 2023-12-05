
function[AD,x,b,info] = DL4IP_ProblemSetup(options)

% set seed for reproducibility
rng(options.data.seed)

% ----------------------------------------------------------------------- %
% load dictionary
out       = load(options.dictionary.filename);

if contains(options.dictionary.filename,'ScheppLogan')
    out.D = out.bestFit.D;
    out.options.data.patchSize = [out.p, out.q];
end
D         = out.D;
patchSize = out.options.data.patchSize;

% ----------------------------------------------------------------------- %
% load data to reconstruct
switch options.data.name
    case 'moon'
        img = imread('moon.tif');
        img = prepareImage(img,options);
        X   = double(img);

    case 'GirlWithPearlEarring'
        % https://en.wikipedia.org/wiki/Girl_with_a_Pearl_Earring#/media/File:1665_Girl_with_a_Pearl_Earring.jpg
        img = imread('GirlWithPearlEarring.jpeg');
        img = prepareImage(img,options);
        X   = double(img);

    case 'peppers'
        img = imread('peppers.png');
        img = prepareImage(img,options);
        X   = double(img);

    case 'tomography'
        % requires AIRToolsII
        X = phantomgallery('shepplogan',options.data.size(1));

    otherwise 
        error("Choose data: {'moon', 'GirlWithPearlEarring', 'peppers', 'tomography'}")
end

% check for incompatible sizes
d1 = size(X,1);
d2 = size(X,2);
if mod(d1,patchSize(1)) ~= 0, d1 = patchSize(1) * floor(d1 / patchSize(1)); end
if mod(d2,patchSize(2)) ~= 0, d2 = patchSize(2) * floor(d2 / patchSize(1));end
X = imresize(X,[d1,d2]);

% ----------------------------------------------------------------------- %
% choose A operator
switch options.AOperator.name
    case 'deblur'
        [A,b,x,info] = deblurringProblem(X,options.AOperator.params{:});

    case 'indicator'
        [A,b,x,info] = deadPixelProblem(X,options.AOperator.params{:});

    case 'superresolution'
        [A,b,x,info] = superResolutionProblem(X,options.AOperator.params{:});
        info.sizeX   = size(X);
        info.sizeB   = [10,10];
        
    case 'denoise'
        [A,b,x,info] = denoiseProblem(X,options.AOperator.params{:});

    case 'tomography'
        % requires AIRToolsII
        [A,b,x,info] = tomographyProblem(X,options.AOperator.params{:});
    
    otherwise
        error("Choose AOperator: {'deblur', 'indicator', 'superresolution', 'denoise', 'tomography'}")
end


% ----------------------------------------------------------------------- %
% choose L operator

switch options.LOperator.name
    case 'Tikhonov'
        L = speye(size(A,2));

    case 'finite difference'
        L = dOperator('finite difference', size(X));

    case 'patch smoother'
        L = dOperatorPatchSmoother(patchSize,size(X));

    case 'patch smoother 2'
        L = dOperatorPatchSmoother2(patchSize,size(X));
    
    case 'none'
        L = [];

    otherwise
        error("Choose LOperator: {'none', 'Tikhonov', 'finite difference', 'patch smoother'}")
end

% ----------------------------------------------------------------------- %
% create A + Dictionary operator
AD = dOperatorDictionaryRepresentation(A,D,patchSize,size(X),options.LOperator.lambda,L);

% add additional info fields (including image resizing)
info.patchSize = patchSize;
info.imgSize   = size(X);

end

%% helper functions
function[img] = prepareImage(img,options)
    if ~isempty(options.data.size)
        img = imresize(img,options.data.size);
    end

    if ~options.data.color && size(img,3) == 3
        img = rgb2gray(img);
    end
end


