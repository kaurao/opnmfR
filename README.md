**Warning** *opnmfR* is still under development and some parts are not completely tested.

# Installation

In R terminal

```
install.packages("remotes")
remotes::install_github("kaurao/opnmfR")
```

**Note** [Microsoft R Open](https://mran.microsoft.com/open) comes with MKL and can substantially improve runtime performance of Rcpp functions by using threads.

# Use it

`data("iris")`

`nn = opnmfR_ranksel_perm(data.matrix(iris[,1:4]), 1:4, W0="nndsvd")`


## Using gpuR (non tested)

Install gpuR following the instrctions here: 

Note that some of the instructions are opudated, e.g. nvidia toolkit for ubuntu14.04. Go to the nvidia page and follow the instructions there.

I needed to do the following on ubuntu 16.04:

`sudo apt-get install opencl-headers ocl-icd-opencl-dev`

Get the correct `.deb` from https://developer.nvidia.com/cuda-downloads and follow the associated instructions.

Run `sudo R` and then `install.packages('gpuR')`



And then run R in sudo mode







