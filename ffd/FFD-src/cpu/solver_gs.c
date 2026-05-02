///////////////////////////////////////////////////////////////////////////////
///
/// \file   solver_gs.c
///
/// \brief  Gauss-Seidel solvers
///
/// \author Mingang Jin, Qingyan Chen
///         Purdue University
///         Jin55@purdue.edu, YanChen@purdue.edu
///         Wangda Zuo
///         University of Miami
///         W.Zuo@miami.edu
///         Wei Tian
///         University of Miami, Schneider Electric
///         w.tian@umiami.edu, Wei.Tian@Schneider-Electric.com
///
/// \date   6/15/2017
///
///////////////////////////////////////////////////////////////////////////////

#include "solver_gs.h"

///////////////////////////////////////////////////////////////////////////////
/// Gauss-Seidel scheme
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param Type Type of variable
///\param x Pointer to variable
///
///\return Residual
///////////////////////////////////////////////////////////////////////////////
int GS_itr(PARA_DATA *para, REAL **var, REAL *x, REAL *flag, REAL residual_min) {
  REAL *as = var[AS], *aw = var[AW], *ae = var[AE], *an = var[AN];
  REAL *ap = var[AP], *af = var[AF], *ab = var[AB], *b = var[B];
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
  int i, j, k;
  REAL SOR = 1.0;
  REAL residual = 0.0;
  int max_iter = 1000; // set the maximum iteration number
  int iter = 0; // iteration counter

  /****************************************************************************
  | OLD: Solve the space using G-S sovler for 5 * 6 = 30 times
  | NEW: Solve the space using G-S sovler until residual is as specified
  ****************************************************************************/
  
  while (1) {
      iter++;

	  // for (it = 0; it < num_swipe; it++) {   // Previous implementation: fixed number of sweeps
          /*-------------------------------------------------------------------------
          | Solve in X in forward direction
          -------------------------------------------------------------------------*/
          for (i = 1; i <= imax; i++)
              for (j = 1; j <= jmax; j++)
                  for (k = 1; k <= kmax; k++) {
                      if (flag[IX(i, j, k)] >= 0) continue;

                      x[IX(i, j, k)] = (1 - SOR) * x[IX(i, j, k)] + SOR * ((ae[IX(i, j, k)] * x[IX(i + 1, j, k)]
                          + aw[IX(i, j, k)] * x[IX(i - 1, j, k)]
                          + an[IX(i, j, k)] * x[IX(i, j + 1, k)]
                          + as[IX(i, j, k)] * x[IX(i, j - 1, k)]
                          + af[IX(i, j, k)] * x[IX(i, j, k + 1)]
                          + ab[IX(i, j, k)] * x[IX(i, j, k - 1)]
                          + b[IX(i, j, k)]) / ap[IX(i, j, k)]);
                      /*if (isnan(x[IX(40, 40, 41)] - x[IX(40, 40, 40)])) {
                          sprintf(msg, "CHECKING_GS %f %f", x[IX(40, 40, 41)] - x[IX(40, 40, 40)], b[IX(i, j, k)]);
                          ffd_log(msg, FFD_NORMAL);
                      }*/
                      /*if (ae[IX(i, j, k)] > 0.001) {
                          sprintf(msg, "CHECKING_GS %e %e %e %e %e %e %e", b[IX(i, j, k)], ae[IX(i, j, k)] * x[IX(i + 1, j, k)], ab[IX(i, j, k)] * x[IX(i, j, k - 1)], aw[IX(i, j, k)] * x[IX(i - 1, j, k)], an[IX(i, j, k)] * x[IX(i, j + 1, k)], as[IX(i, j, k)] * x[IX(i, j - 1, k)], af[IX(i, j, k)] * x[IX(i, j, k + 1)]);
                          ffd_log(msg, FFD_NORMAL);
                      }*/
                      /*if isnan(x[IX(i, j, k)]) {
                          sprintf(msg, "DIVERGENCE ERROR");
                          ffd_log(msg, FFD_ERROR);
                      }*/
                  }
          /*-------------------------------------------------------------------------
          | Solve in X in backward direction
          -------------------------------------------------------------------------*/
          for (i = imax; i >= 1; i--)
              for (j = 1; j <= jmax; j++)
                  for (k = 1; k <= kmax; k++) {
                      if (flag[IX(i, j, k)] >= 0) continue;

                      x[IX(i, j, k)] = (1 - SOR) * x[IX(i, j, k)] + SOR * ((ae[IX(i, j, k)] * x[IX(i + 1, j, k)]
                          + aw[IX(i, j, k)] * x[IX(i - 1, j, k)]
                          + an[IX(i, j, k)] * x[IX(i, j + 1, k)]
                          + as[IX(i, j, k)] * x[IX(i, j - 1, k)]
                          + af[IX(i, j, k)] * x[IX(i, j, k + 1)]
                          + ab[IX(i, j, k)] * x[IX(i, j, k - 1)]
                          + b[IX(i, j, k)]) / ap[IX(i, j, k)]);
                      /*if isnan(x[IX(i, j, k)]) {
                          sprintf(msg, "DIVERGENCE ERROR");
                          ffd_log(msg, FFD_ERROR);
                      }*/
                  }
          ///*-------------------------------------------------------------------------
          //| Solve in Y in forward direction
          //-------------------------------------------------------------------------*/
          for (j = 1; j <= jmax; j++)
              for (i = 1; i <= imax; i++)
                  for (k = 1; k <= kmax; k++) {
                      if (flag[IX(i, j, k)] >= 0) continue;

                      x[IX(i, j, k)] = (1 - SOR) * x[IX(i, j, k)] + SOR * ((ae[IX(i, j, k)] * x[IX(i + 1, j, k)]
                          + aw[IX(i, j, k)] * x[IX(i - 1, j, k)]
                          + an[IX(i, j, k)] * x[IX(i, j + 1, k)]
                          + as[IX(i, j, k)] * x[IX(i, j - 1, k)]
                          + af[IX(i, j, k)] * x[IX(i, j, k + 1)]
                          + ab[IX(i, j, k)] * x[IX(i, j, k - 1)]
                          + b[IX(i, j, k)]) / ap[IX(i, j, k)]);
                      /*if isnan(x[IX(i, j, k)]) {
                          sprintf(msg, "DIVERGENCE ERROR");
                      }*/
                  }
          ///*-------------------------------------------------------------------------
          //| Solve in Y in backward direction
          //-------------------------------------------------------------------------*/
          for (j = jmax; j >= 1; j--)
              for (i = 1; i <= imax; i++)
                  for (k = 1; k <= kmax; k++) {
                      if (flag[IX(i, j, k)] >= 0) continue;

                      x[IX(i, j, k)] = (1 - SOR) * x[IX(i, j, k)] + SOR * ((ae[IX(i, j, k)] * x[IX(i + 1, j, k)]
                          + aw[IX(i, j, k)] * x[IX(i - 1, j, k)]
                          + an[IX(i, j, k)] * x[IX(i, j + 1, k)]
                          + as[IX(i, j, k)] * x[IX(i, j - 1, k)]
                          + af[IX(i, j, k)] * x[IX(i, j, k + 1)]
                          + ab[IX(i, j, k)] * x[IX(i, j, k - 1)]
                          + b[IX(i, j, k)]) / ap[IX(i, j, k)]);
                      /*if isnan(x[IX(i, j, k)]) {
                          sprintf(msg, "DIVERGENCE ERROR");
                          ffd_log(msg, FFD_ERROR);
                      }*/
                  }
          ///*-------------------------------------------------------------------------
          //| Solve in Z in forward direction
          //-------------------------------------------------------------------------*/
          for (k = 1; k <= kmax; k++)
              for (i = 1; i <= imax; i++)
                  for (j = 1; j <= jmax; j++) {
                      if (flag[IX(i, j, k)] >= 0) continue;

                      x[IX(i, j, k)] = (1 - SOR) * x[IX(i, j, k)] + SOR * ((ae[IX(i, j, k)] * x[IX(i + 1, j, k)]
                          + aw[IX(i, j, k)] * x[IX(i - 1, j, k)]
                          + an[IX(i, j, k)] * x[IX(i, j + 1, k)]
                          + as[IX(i, j, k)] * x[IX(i, j - 1, k)]
                          + af[IX(i, j, k)] * x[IX(i, j, k + 1)]
                          + ab[IX(i, j, k)] * x[IX(i, j, k - 1)]
                          + b[IX(i, j, k)]) / ap[IX(i, j, k)]);
                      /*if isnan(x[IX(i, j, k)]) {
                          sprintf(msg, "DIVERGENCE ERROR");
                          ffd_log(msg, FFD_ERROR);
                      }*/
                  }
          ///*-------------------------------------------------------------------------
          //| Solve in Z in backward direction
          //-------------------------------------------------------------------------*/
          for (k = kmax; k >= 1; k--)
              for (i = 1; i <= imax; i++)
                  for (j = 1; j <= jmax; j++) {
                      if (flag[IX(i, j, k)] >= 0) continue;

                      x[IX(i, j, k)] = (1 - SOR) * x[IX(i, j, k)] + SOR * ((ae[IX(i, j, k)] * x[IX(i + 1, j, k)]
                          + aw[IX(i, j, k)] * x[IX(i - 1, j, k)]
                          + an[IX(i, j, k)] * x[IX(i, j + 1, k)]
                          + as[IX(i, j, k)] * x[IX(i, j - 1, k)]
                          + af[IX(i, j, k)] * x[IX(i, j, k + 1)]
                          + ab[IX(i, j, k)] * x[IX(i, j, k - 1)]
                          + b[IX(i, j, k)]) / ap[IX(i, j, k)]);
                      /*if isnan(x[IX(i, j, k)]) {
                          sprintf(msg, "DIVERGENCE ERROR");
                          ffd_log(msg, FFD_ERROR);
                      }*/
                  }

      //} END OF PREVIOUS IMPLEMENTATION OF FIXED NUMBER OF SWIPES

      //QFLUXFLUX
      residual = check_residual(para, var, x, flag);
	  //sprintf(msg, "EWAN: Residual: %f", residual); // DEBUG
      //ffd_log(msg, FFD_NORMAL); // DEBUG
      if (residual < residual_min || iter == max_iter) break;
  }

  residual = check_residual(para, var, x, flag);
  sprintf(msg, "Number of Iterations: %d", iter);
  ffd_log(msg, FFD_NORMAL);
 
 /* sprintf(msg, "Residual in the GS solver is: %e", residual);
  ffd_log(msg, FFD_NORMAL);*/
  //printf("residual in the solver is %f\n", check_residual(para, var, x));
  return 0;
} // End of GS_itr()

  
///////////////////////////////////////////////////////////////////////////////
/// Gauss-Seidel solver
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param flag Pointer to the cell property flag
///\param x Pointer to variable
///
///\return Residual
///////////////////////////////////////////////////////////////////////////////
int Gauss_Seidel(PARA_DATA *para, REAL **var,  REAL *x, REAL *flag, REAL residual) {
  GS_itr(para, var,  x, flag, residual);
  return 0;
} // End of Gauss-Seidel( )


///////////////////////////////////////////////////////////////////////////////
/// Jacobi Scheme for pressure
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param Type Type of variable
///\param x Pointer to variable
///
///\return Residual
///////////////////////////////////////////////////////////////////////////////
int Jacobi_iter(PARA_DATA *para, REAL **var, REAL *x,REAL *flag, int num_swipe) {
  REAL *as = var[AS], *aw = var[AW], *ae = var[AE], *an = var[AN];
  REAL *ap = var[AP], *af = var[AF], *ab = var[AB], *b = var[B];
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
  int i, j, k, it;
  REAL *tmp = var[TMP4];
  REAL *flagp=var[FLAGP];

  /****************************************************************************
  | Solve the space using Jacobi sovler for num_swipe * 6 = 30 times
  ****************************************************************************/
  //while (residual > 1e-6) {
  for (it = 0; it<num_swipe*6; it++) {
    /*-------------------------------------------------------------------------
    | Solve in X(1->imax), Y(1->jmax), Z(1->kmax)
    -------------------------------------------------------------------------*/
    for (i = 1; i <= imax; i++)
      for (j = 1; j <= jmax; j++)
        for (k = 1; k <= kmax; k++) {
          if (flag[IX(i, j, k)] >= 0 ) continue;


          tmp[IX(i, j, k)] = (ae[IX(i, j, k)] * x[IX(i + 1, j, k)]
            + aw[IX(i, j, k)] * x[IX(i - 1, j, k)]
            + an[IX(i, j, k)] * x[IX(i, j + 1, k)]
            + as[IX(i, j, k)] * x[IX(i, j - 1, k)]
            + af[IX(i, j, k)] * x[IX(i, j, k + 1)]
            + ab[IX(i, j, k)] * x[IX(i, j, k - 1)]
            + b[IX(i, j, k)]) / ap[IX(i, j, k)];
        }

    for (i = 1; i <= imax; i++)
      for (j = 1; j <= jmax; j++)
        for (k = 1; k <= kmax; k++) {
          if (flag[IX(i, j, k)] >= 0 ) continue;

          x[IX(i, j, k)] = tmp[IX(i, j, k)];

        }

  }

  return 0;
} // End of Jacobi_P()


///////////////////////////////////////////////////////////////////////////////
/// Jacobi solver
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param flag Pointer to the cell property flag
///\param x Pointer to variable
///
///\return Residual
///////////////////////////////////////////////////////////////////////////////
int Jacobi(PARA_DATA *para, REAL **var, REAL *flag, REAL *x, int num_swipe) {
  Jacobi_iter(para, var, x, flag, num_swipe);
  return 0;

} // End of Jacobi( )
