#!/bin/bash

# document R code
R -f roxygen.R
R CMD Rd2pdf --force --no-preview -o manual.pdf ./pkg

# document the C API
if ! [ -x "$(command -v doxygen)" ]; then
  echo 'Warning: Doxygen is not installed. Exiting' >&2
  exit 0
fi


basedir=`pwd`
cd pkg/inst/include
doxygen Doxyfile
cd $basedir
cd pkg/inst/doc/latex 
make
cd ..
mv latex/refman.pdf ./stringdist_api.pdf
rm -rf latex
cd $basedir


