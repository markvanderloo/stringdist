#!/bin/bash

# perform build and R CMD check with undefined behaviour sanitizer switched on.

## Thanks to Dirk for posting this command here:
## http://dirk.eddelbuettel.com/blog/2015/01/18/#ubsan-clang-container

cd output

R CMD build ../pkg 

docker run --rm -ti -v $(pwd):/mnt rocker/r-devel-ubsan-clang check.r --setwd /mnt --install-deps stringdist_*.tar.gz

cd ..

