#!/bin/bash

# document R code
R -e "pkgload::load_all('pkg'); roxygen2::roxygenize('pkg')"
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
cd pkg/vignettes/latex 
make
cd ..
mv latex/refman.pdf ./stringdist_api.pdf
rm -rf latex
cd $basedir


