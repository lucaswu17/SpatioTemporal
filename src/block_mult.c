#include "R.h"
#include "Rinternals.h"
#include <R_ext/Lapack.h>

//multiplication of block diagonal matrix with another matrix
//inputs are:
//matrix size (BLOCK matrix)
//matrix size (multiplication target)
//number of blocks
//vector with the size of each block
//the block-matrix it self
//the matrix to multiply with
// BLOCK*X
SEXP block_mult(SEXP n_totR, SEXP n_xR, SEXP n_blocksR, 
		SEXP block_sizesR, SEXP matR, SEXP XR){

  int i, k, l, j, n_cum;
  int n_blocks, *block_sizes, n_tot, n_x;
  double *mat, *X, *result;
  SEXP resultR;

  n_tot = *INTEGER(n_totR);
  n_x = *INTEGER(n_xR);
  n_blocks = *INTEGER(n_blocksR);
  block_sizes = INTEGER(block_sizesR);
  mat = REAL(matR);
  X = REAL(XR);
  PROTECT(resultR = allocMatrix(REALSXP, n_tot, n_x));
  result = REAL(resultR);

  //set results to zero
  memset(result, 0, n_tot*n_x*sizeof(double));
    
  for(j=0; j<n_x; ++j){
    n_cum = 0;
    for(i=0; i < n_blocks; ++i){
      for(k=0; k < block_sizes[i]; ++k)
	for(l=0; l < block_sizes[i]; ++l)
	  result[n_cum+k+j*n_tot] += mat[n_cum+k + (n_cum+l)*n_tot] * 
	    X[n_cum+l+j*n_tot];
      n_cum += block_sizes[i];
    }
  }

  UNPROTECT(1); /* resultR */
  return resultR;
}
