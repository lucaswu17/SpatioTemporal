#include "R.h"
#include "Rinternals.h"
#include <R_ext/Lapack.h>

//declaration of external FORTRAN function
extern void F77_CALL(dpodi)(double*, int*, int*, double*, int*);

//calculate a mtrix inverse from a cholesky factor (requires a previous call
//to make_chol_block)
//inputs are:
//matrix size
//number of blocks
//vector with the size of each block
//maximum size of any block
//the matrix it self.
SEXP inv_chol_block(SEXP n_totR, SEXP n_blocksR, SEXP block_sizesR, SEXP max_sizeR, SEXP matR){

  int i, j, k, n_this, n_cum;
  int job=1; //compute inverse only (no determinant)
  int n_blocks, *block_sizes, max_size, n_tot;
  double *mat, *result, *block, det;
  SEXP resultR;

  n_tot = *INTEGER(n_totR);
  n_blocks = *INTEGER(n_blocksR);
  block_sizes = INTEGER(block_sizesR);
  max_size = *INTEGER(max_sizeR);
  mat = REAL(matR);
  PROTECT(resultR = allocMatrix(REALSXP, n_tot, n_tot));
  result = REAL(resultR);

  if(n_blocks==1){ //don't need to copy matrix if only one block.
    memcpy(result, mat, n_tot*n_tot*sizeof(double));
	
    F77_CALL(dpodi)(result, &n_tot, &n_tot, &det, &job);
    for (i=1; i < n_tot; ++i){
      for(j=0; j < i; ++j){
	result[i + j*n_tot] = result[j + i*n_tot];
      }
    }
  }else{
    block = malloc(max_size*max_size * sizeof (double));
    memset(result, 0, n_tot*n_tot*sizeof(double));
	
    n_cum = 0;
    for(i=0; i < n_blocks; ++i){
      n_this = block_sizes[i];
      for (k=0; k < n_this; ++k){
	for(j=0; j <=k; ++j){
	  block[j+n_this*k] = mat[j+n_cum + (k+n_cum)*n_tot];
	}
      }
      F77_CALL(dpodi)(block, &n_this, &n_this, &det, &job);
      for (k=0; k < n_this; ++k){
	for(j=0; j <= k; ++j){
	  result[j+n_cum + (k+n_cum)*n_tot] = block[j+n_this*k]; 
	  result[k+n_cum + (j+n_cum)*n_tot] = block[j+n_this*k]; 
	}
      }
      n_cum += n_this;
    }      
    free(block);
  }
  
  UNPROTECT(1); /* resultR */
  return resultR;
}
