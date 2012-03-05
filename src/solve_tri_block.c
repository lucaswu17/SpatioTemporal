#include "R.h"
#include "Rinternals.h"
#include <R_ext/Lapack.h>

//declaration of external FORTRAN function
extern void F77_CALL(dtrsl)(double*, int*, int*, double*, int*, int*);

//solves the system T*X=B where T is block upper triagular (output from
//make_chol_block)
//inputs is total matrix side, number of blocks, vector with the size of each
//block, maximum size of any block, number of system to solve (dim(B)[2]),
//indicator wether or not to transpose T, the triangular matrix, and the
//LHS (B)of the equation system.
SEXP solve_tri_block(SEXP n_totR, SEXP n_blocksR, SEXP block_sizesR, SEXP max_sizeR, SEXP n_xR, SEXP transposeR, SEXP matR, SEXP xR){

  int i, j, k, n_this, n_cum, n_x;
  int info, job;
  int n_blocks, *block_sizes, max_size, n_tot;
  double *x, *mat, *result, *block;
  SEXP resultR;

  n_tot = *INTEGER(n_totR);
  n_blocks = *INTEGER(n_blocksR);
  block_sizes = INTEGER(block_sizesR);
  max_size = *INTEGER(max_sizeR);
  n_x = *INTEGER(n_xR);
  //job defines T*X=B to solve. T is upper triangular (=1)
  //and possibly transpose (+10)
  job = 10*(*INTEGER(transposeR)!=0)+1; 
  mat = REAL(matR);
  x = REAL(xR);
  PROTECT(resultR = allocMatrix(REALSXP, n_tot, n_x));
  result = REAL(resultR);

  memcpy(result, x, n_tot*n_x*sizeof(double));

  if(n_blocks==1){ //don't need to copy matrix if only one block.
    for(i=0; i<n_x; ++i){
      F77_CALL(dtrsl)(mat, &n_tot, &n_tot, result+i*n_tot, &job, &info);
    }
  }else{
    block = malloc(max_size*max_size * sizeof (double));
    n_cum = 0;
    for(k=0; k < n_blocks; ++k){
      n_this = block_sizes[k];
      for (i=0; i < n_this; ++i){
	for(j=0; j <=i; ++j){
	  block[j+n_this*i] = mat[j+n_cum + (i+n_cum)*n_tot];
	}
      }
      for(i=0; i<n_x; ++i){
	F77_CALL(dtrsl)(block, &n_this, &n_this, result+n_cum+i*n_tot, &job, &info);
      }
      n_cum += n_this;
    }
    free(block);
  }
  UNPROTECT(1); /* resultR */
  return resultR;
}
