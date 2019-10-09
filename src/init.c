#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>


void DAG_int_var(double *X, double *A, double *B, double *C, double *U, int *mm, int *MM, int *nn, 
    double *Lambda, double *Lambda2, double *tautau, int *A_NZ, int *NZ, double *sigma, 
    double *sigmaX, double *Sig, double *tol, double *obj, double *XTX, double *XTX_inv,
    double *rhorho, double *rhorho2, double *rhorho3, int *maxIter);   

void DAG_int(double *X, double *A, double *B, int *mm, int *MM, int *nn, double *Lambda, double *Lambda2,
       double *tautau, int *A_NZ, int *NZ, double *sigma, double *sigmaX, double *Sig, double *tol,
       double *obj, double *XTX, double *XTX_inv, double *X2, double *rhorho, double *rhorho2, int *maxIter);        

void DAG_obs(double *X, double *A, int *mm, int *nn, double *Lambda,
       double *tautau, int *A_NZ, int *NZ, double *sigma, double *tol, double *obj,
       double *XTX, double *XTX_inv, double *rhorho, int *maxIter);


static R_NativePrimitiveArgType DAG_int_var_types[24] =
  {REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP,
  	REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP,
   	REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};

static R_NativePrimitiveArgType DAG_int_types[22] =
  {REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
    REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};

static R_NativePrimitiveArgType DAG_obs_types[15] =
  {REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP,
  	REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};


static R_CMethodDef cMethods[] = {
  {"DAG_int_var", (DL_FUNC) &DAG_int_var, 24, DAG_int_var_types},
  {"DAG_int", (DL_FUNC) &DAG_int, 22, DAG_int_types},
  {"DAG_obs", (DL_FUNC) &DAG_obs, 15, DAG_obs_types},
  {NULL, NULL, 0, NULL}
};
 

void attribute_visible R_init_intdag(DllInfo *info) {
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
