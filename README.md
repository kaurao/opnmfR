# Installation

In the system terminal

`cd ~`

`git clone https://jugit.fz-juelich.de/inm7/aml/opnmfr`


In R terminal

`library(devtools)`

`devtools::install('~/opnmfr')`


# Use it

`data("iris")`

`nn = opnmfR_ranksel_perm(data.matrix(iris[,1:4]), 1:4, W0="nndsvd")`


## Using gpuR

Install gpuR following the instrctions here: 

Note that some of the instructions are opudated, e.g. nvidia toolkit for ubuntu14.04. Go to the nvidia page and follow the instructions there.

I needed to do the following on ubuntu 16.04:

`sudo apt-get install opencl-headers ocl-icd-opencl-dev`

Get the corect `.deb` from here https://developer.nvidia.com/cuda-downloads and follow the associated instructions.

Run `sudo R` and then `install.packages('gpuR')`



And then run R in sudo mode







