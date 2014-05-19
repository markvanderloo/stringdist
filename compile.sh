#!/bin/bash

cd pkg/src
rm --verbose *.o *.so
gcc -std=gnu99 -I/usr/share/R/include -DNDEBUG -fpic  -O0 -Wall -pipe  -g  -c *.c 
gcc -std=gnu99 -shared -o stringdist.so *.o -L/usr/lib/R/lib -lR
cd ../../

