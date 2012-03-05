#include "R.h"
#include "Rinternals.h"
#include <R_ext/Lapack.h>
#include <string.h>

//declaration of external FORTRAN function
extern void F77_CALL(dpofa)(double*, int*, int*, int*);

//cholesky inverse of a block diagonal matrix using the diagonal structure.
//inputs are:
//matrix size
//number of blocks
//vector with the size of each block
//maximum size of any block
//the matrix it self.
SEXP make_chol_block(SEXP n_totR, SEXP n_blocksR, SEXP block_sizesR,
		     SEXP max_sizeR, SEXP matR){

  int i, j, k, n_this, n_cum;
  int info,pos_def=1;
  int n_blocks, *block_sizes, max_size, n_tot;
  double *mat, *result, *block;
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
    
    F77_CALL(dpofa)(result, &n_tot, &n_tot, &info);
    //info info!=0 then the matrix is not positive definite
    pos_def = (info==0);
    //set lower block to zero
    for (i=0; i < n_tot; ++i){
      for(j=i+1; j < n_tot; ++j){
	result[j + i*n_tot] = 0;
      }
    }
  }else{
    block = malloc(max_size*max_size * sizeof (double));
    memset(result, 0, n_tot*n_tot*sizeof(double));
    
    n_cum = 0;
    for(i=0; i < n_blocks; i++){
      n_this = block_sizes[i];
      for(k=0; k < n_this; ++k){
	for(j=0; j <=k; ++j){
	  block[j+n_this*k] = mat[j+n_cum + (k+n_cum)*n_tot];
	}
      }
      F77_CALL(dpofa)(block, &n_this, &n_this, &info);
      //info info!=0 then the matrix is not positive definite
      pos_def = (info==0);
      for (k=0; k < n_this; ++k){
	for(j=0; j <= k; ++j){
	  result[j+n_cum + (k+n_cum)*n_tot] = block[j+n_this*k]; 
	}
      }    
      n_cum += n_this;
    }
    free(block);
  }
  if(!pos_def) //not positive definite, set top left value to -1
    result[0 + 0*n_tot] = -1;

  UNPROTECT(1); /* resultR */
  return resultR;
}
