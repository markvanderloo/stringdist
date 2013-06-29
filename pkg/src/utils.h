
#ifndef sd_utils_h
#define sd_utils_h

/* integer recycling macro ( equals X % Y) */
#define RECYCLE(X,Y) ( (X) == (Y) ? 0 : (X) )


double min3(double, double, double);

double min2(double, double);

unsigned int max_length(SEXP);


#endif
