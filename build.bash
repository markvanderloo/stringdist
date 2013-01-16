#!/bin/bash

R=R
if [ ${#} -gt 0 ]; then
    if [ "$1" = "-dev" ]; then
        R=Rdev  
    fi
fi

echo "######## Removing building information..."
rm -rf output

echo "######## Copying DESCRIPTION and NAMESPACE to pkg directory..."
cp build/DESCRIPTION pkg
cp build/NAMESPACE pkg

echo "######## Generate documentation..."
$R -q -f roxygen.R

echo "######## Building package in output..."
mkdir output
cd output
$R CMD build ../pkg
echo "######## Testing package..."
for x in *.tar.gz 
do 
    $R CMD check --as-cran $x
done

echo "**BUILT USING $R"
$R --version

