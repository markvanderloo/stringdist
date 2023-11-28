#!/bin/bash

cd pkg/src
rm --verbose *.o *.so
gcc-10 -std=gnu99 -I/usr/share/R/include -DNDEBUG -fpic -fopenmp -O4 -Wall -pipe  -g  -c -Wstrict-prototypes -Wformat *.c 
#gcc -std=gnu99 -I/usr/share/R/include -DNDEBUG -fpic -O2 -Wall -pipe  -g  -c *.c 
gcc -std=gnu99 -shared -o stringdist.so *.o -L/usr/lib/R/lib -lR
cd ../../

