# DL4IP
Dictionary Learning for Inverse Problems

## Installation and Dependencies

To gain access to the code, download or clone from Github using
```console
git clone https://github.com/elizabethnewman/DL4IP.git
```
This will install the directory ```DL4IP``` on your machine. Open MATLAB, change to the ```DL4IP``` directory, and run the setup file:
```console
cd DL4IP
DL4IP_PathSetup
```
This path setup file will automatically clone the dependent repositories. 
* [IRTools](https://github.com/jnagy1/IRtools): orginal MRNSD algorithm and 
* [AIRTools](https://github.com/jakobsj/AIRToolsII): problem setup for certain examples
* [matlab2tikz](https://github.com/matlab2tikz/matlab2tikz): way to generate nice plots




## Repository Structure

The repository has the following structure:
* [code](https://github.com/elizabethnewman/DL4IP/tree/main/code): main code for solving inverse problems using sparse dictionary representations
* [dictionaries](https://github.com/elizabethnewman/DL4IP/tree/main/dictionaries): pre-formed dictionaries from [built-in MATLAB Flowers dataset](https://www.mathworks.com/help/deeplearning/ug/data-sets-for-deep-learning.html) and [high-resolution image of the earth](https://www.cnet.com/science/stunning-high-resolution-photo-shows-earths-many-hues/).
* [examples](https://github.com/elizabethnewman/DL4IP/tree/main/examples): scripts to generate the examples in our paper.


## Citation

```
@misc{newman2023image,
      title={Image reconstructions using sparse dictionary representations and implicit, non-negative mappings}, 
      author={Elizabeth Newman and Jack Michael Solomon and Matthias Chung},
      year={2023},
      eprint={2312.03180},
      archivePrefix={arXiv},
      primaryClass={math.NA}
}
```



