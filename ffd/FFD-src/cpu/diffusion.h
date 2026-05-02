///////////////////////////////////////////////////////////////////////////////
///
/// \file   diffusion.c
///
/// \brief  Calculate the diffusion equation
///
/// \author Wangda Zuo
///         University of Miami
///         W.Zuo@miami.edu
///         Mingang Jin, Qingyan Chen
///         Purdue University
///         Jin55@purdue.edu, YanChen@purdue.edu
///         Wei Tian
///         University of Miami, Schneider Electric
///         w.tian@umiami.edu, Wei.Tian@Schneider-Electric.com
///
/// \date   6/15/2017
///
/// This file provides functions that are used for calculating the diffusion
/// equations.
///
///////////////////////////////////////////////////////////////////////////////
#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include "data_structure.h"

#include "boundary.h"

#include "solver.h"

#include "utility.h"

#include "chen_zero_equ_model.h"


///////////////////////////////////////////////////////////////////////////////
/// Entrance of calculating diffusion equation
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type Type of variable
///\param index Index of trace substance or species
///\param psi Pointer to the variable at current time step
///\param psi0 Pointer to the variable at previous time step
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int diffusion(PARA_DATA *para, REAL **var, int var_type, int index,
               REAL *psi, REAL *psi0, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Calculate coefficients for diffusion equation solver
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param psi Pointer to the variable at current time step
///\param psi0 Pointer to the variable at previous time step
///\param var_type Type of variable
///\param index Index of trace substance or species
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int coef_diff(PARA_DATA *para, REAL **var, REAL *psi, REAL *psi0,
               int var_type, int index, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Calculate source term in the diffusion equation
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type Type of variable
///\param index Index of trace substances or species
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int source_diff(PARA_DATA *para, REAL **var, int var_type, int index);

#endif
