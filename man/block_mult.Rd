\name{block_mult}
\encoding{latin1}
\Rdversion{1.1}
\alias{block_mult}
\alias{calc_F_part_X}
\alias{calc_tF_mat_F}
\alias{calc_tF_times_mat}
\alias{dot_prod}
\alias{inv_chol_block}
\alias{make_chol_block}
\alias{make_sigma_b}
\alias{make_sigma_b_full}
\alias{make_sigma_nu}
\alias{make_sigma_nu_cross_cov}
\alias{norm2_c}
\alias{solve_tri_block}
\alias{sum_log_diag}

\title{
  Internal C-functions
}
\description{
  Handles of the internal C-functions. \strong{DO NOT USE}, call the
  wrappers instead.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom} and Adam Szpiro
}
\section{Warning}{
  These are internal C-functions, incorrect usage \strong{will crash}
  R, use the provided wrappers instead.
}
\seealso{
  Functions that provide R interface/wrappers for the C-code: \cr
  \code{\link{dot.prod}}, \code{\link{norm2}}, \code{\link{sumLogDiag}},
  \code{\link{makeCholBlock}}, \code{\link{invCholBlock}},
  \code{\link{solveTriBlock}}, \code{\link{make.sigma.B}},
  \code{\link{make.sigma.B.full}}, \code{\link{make.sigma.nu}}, \cr
  \code{\link{make.sigma.nu.cross.cov}},
  \code{\link{calc.tF.times.mat}}, \code{\link{calc.tF.mat.F}}, and \cr
  \code{\link{calc.F.times.X}}.
}
