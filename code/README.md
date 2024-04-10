## Main Setup

Our goal is to solve inverse problems using sparse dictionary representations. 


## Main Functions

The main functions for running this code include the prefix ```DL4IP_```.  The functions are as follows:
* ```DL4IP_OptionsDefault```: choose options for running experiment
  * ```options.data```: parameters for setting up the data
  * ```options.AOperator```: choose options for problem forward operator ('deblur', 'indicator', 'superresolution', 'tomography'). The dictionary representation will be added automatically.
  * ```options.LOpertor```: choose options for the regularization operator on the image reconstruction ('none', 'Tikhonov', 'finite difference', 'patch smoother')
  * ```options.optimizer```: choose optimizer ('MRNSDSparsity', 'GDNN')
  * 
* ```DL4IP_DataSetup```: choose the original image to reconstruct. 
