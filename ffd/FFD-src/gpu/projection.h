///////////////////////////////////////////////////////////////////////////////
///
/// \file   projection.c
///
/// \brief  Solver for projection step
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
#ifndef _PROJECTION_H
#define _PROJECTION_H

#include "data_structure.h"

#include "solver_gs.h"

#include "utility.h"

#include "boundary.h"

///////////////////////////////////////////////////////////////////////////////
/// Project the velocity
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int project(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Check the mass imbalance after projjection
/// This usually indicates that the energy balance after advection could be problematic.
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int check_mass_imbalance(PARA_DATA *para, REAL **var);

#endif
