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



