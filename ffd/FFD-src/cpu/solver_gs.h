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

#ifndef _SOLVER_GS_H
#define _SOLVER_GS_H

#include "data_structure.h"

#include "boundary.h"

#include "utility.h"


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
int GS_itr(PARA_DATA *para, REAL **var, REAL *x, REAL *flag, int num_swipe);

///////////////////////////////////////////////////////////////////////////////
/// Gauss-Seidel solver for terms other than pressure
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param flag Pointer to the cell property flag
///\param x Pointer to variable
///
///\return Residual
///////////////////////////////////////////////////////////////////////////////
int Gauss_Seidel(PARA_DATA *para, REAL **var, REAL *x, REAL *flagp, int num_swipe);

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
int Jacobi_iter(PARA_DATA *para, REAL **var, REAL *x, REAL *flag, int num_swipe);

///////////////////////////////////////////////////////////////////////////////
/// Jacobi solver for pressure
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param Type Type of variable
///\param x Pointer to variable
///
///\return Residual
///////////////////////////////////////////////////////////////////////////////
int Jacobi(PARA_DATA *para, REAL **var, REAL *flag, REAL *x, int num_swipe);

#endif
