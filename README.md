# BAS Master Sample
Package for generating sample designs based on a Balanced Acceptance Sampling (BAS) master sample 
for general monitoring. The package is specifically designed for for Western Canada Marine Master Sample
purposes.


# Example:
devtools::install_github("paul-vdb/DFO-master-sample")

library(BASMasterSample)

library(sf)

library(sp)

data(Fed_MPAs_clipped)

smp <- masterSample(Fed_MPAs_clipped, N = 100)

plot(smp)
