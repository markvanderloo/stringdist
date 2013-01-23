// we're comparing 2 or 3 doubles in several places.
static inline double min3(double x, double y, double z){
   if ( x <= y && x <= z ){ 
      return(x);
   } else if ( y <= x && y <= z ) {
      return(y);
   } else {
      return(z);
   }
}

static inline double min2(double x, double y){
   if ( x <= y ){
      return(x);
   } else {
      return(y);
   }
}



