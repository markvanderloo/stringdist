#!/bin/bash

R -f roxygen.R
R CMD Rd2pdf --force --no-preview -o manual.pdf ./pkg


