The naming convention of each example from our paper is the following:
```
ex<PROBLEM><ALGORITHM>
```

For example, to generate the results using gradient descent with a nonnegativity mapping from Figure 5 (image deblurring) in the paper, run the following script in your Matlab console:
```
exDeblurGDNN
```


We offer four different problem types plus additional parameter tuning and and two algorithms.  

**PROBLEM**
* ```Deblur```
* ```Indicator```
* ```Superresolution```
* ```Tomography```
* ```ParameterTuning```

**ALGORITHM**
* ```GDNN```: Gradient Descent with Non-Negativity Mapping
* ```MRNSDSparsity```: Modified Residual Norm Steepest Descent with Sparsity-Promoting Mapping

