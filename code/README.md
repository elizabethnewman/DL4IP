## Main Setup

Our goal is to solve inverse problems using sparse dictionary representations. We solve
$$\min_{\mathbf{x} \ge 0} \|\mathbf{A} \mathbf{x} - \mathbf{b}\|_2^2$$

The scripts to run experiments approximate the following structure:

```
% setup default options
options = DL4IP_OptionsDefault();

% specify problem-dependent options in script

% set seed for reproducibility
rng(options.seed)

% setup problem
[AD,x,b,probInfo] = DL4IP_ProblemSetup(options);

% run algorithms
[coeff,optInfo] = ALGORITHM(AD, [b(:);zeros(size(AD.L,1),1)],varargin,options.ALGORITHM.options);

% compute results
results = computeResults(coeff,optInfo,probInfo,x,AD.D,endTime,'GDNN',options);

```

## Main Functions

The main functions for running this code include the prefix ```DL4IP_```.  The functions are as follows:
* ```DL4IP_ProblemSetup```: setup the problem based on the chosen options. This function does the following:
  * Loads dictionary
  * Loads and preprocesses original image
  * Constructs ```AOperator```and right-hand side ```b``` based on problem (using ```IRTools```, ```AIRTools```, or proprietary problem setup)
  * Constructs regularization operator ```LOperator```)
  * Creates one operator with dictionary representation and regularization using the ```dOperatorDictionaryRepresentation``` class.

* ```DL4IP_OptionsDefault```: choose options for running experiment
  * ```options.data```: parameters for setting up the data
  * ```options.AOperator```: choose options for problem forward operator ('deblur', 'indicator', 'superresolution', 'tomography'). The dictionary representation will be added automatically.
  * ```options.LOpertor```: choose options for the regularization operator on the image reconstruction ('none', 'Tikhonov', 'finite difference', 'patch smoother')
  * ```options.optimizer```: choose optimizer ('MRNSDSparsity', 'GDNN')
  * ```options.dictionary```: choose dictionary to use for reconstruction
  * ```options.seed```: seed for reproducibility (default = 42)


* ```DL4IP_DataSetup```: choose the original image to reconstruct; pairs with ```options.data```, which include
  * ```options.data.name```: image name ('moon', 'GirlWithPearlEarring', 'peppers', 'tomography')
  * ```options.data.size```: size of image (should be divisible by dictionary patch size)
  * ```options.data.color```: Boolean to keep image in color or turn to gray scale
  * ```options.data.seed```: seed for data for reproducibility (e.g., adding noise)

