#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP R_all_int(SEXP);
extern SEXP R_amatch(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_get_qgrams(SEXP, SEXP);
extern SEXP R_lengths(SEXP);
extern SEXP R_lower_tri(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_soundex(SEXP, SEXP);
extern SEXP R_stringdist(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"R_all_int",    (DL_FUNC) &R_all_int,     1},
    {"R_amatch",     (DL_FUNC) &R_amatch,     12},
    {"R_get_qgrams", (DL_FUNC) &R_get_qgrams,  2},
    {"R_lengths",    (DL_FUNC) &R_lengths,     1},
    {"R_lower_tri",  (DL_FUNC) &R_lower_tri,   8},
    {"R_soundex",    (DL_FUNC) &R_soundex,     2},
    {"R_stringdist", (DL_FUNC) &R_stringdist,  9},
    {NULL, NULL, 0}
};

void R_init_stringdist(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
    
    /* used by external packages linking to internal xts code from C */
    R_RegisterCCallable("stringdist","R_all_int",(DL_FUNC) &R_all_int);
    R_RegisterCCallable("stringdist","R_amatch",(DL_FUNC) &R_amatch);
    R_RegisterCCallable("stringdist","R_get_qgrams",(DL_FUNC) &R_get_qgrams);
    R_RegisterCCallable("stringdist","R_lengths",(DL_FUNC) &R_lengths);
    R_RegisterCCallable("stringdist","R_lower_tri",(DL_FUNC) &R_lower_tri);
    R_RegisterCCallable("stringdist","R_soundex",(DL_FUNC) &R_soundex);
    R_RegisterCCallable("stringdist","R_stringdist",(DL_FUNC) &R_stringdist);
}
