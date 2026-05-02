///////////////////////////////////////////////////////////////////////////////
///
/// \file   initialization.c
///
/// \brief  Set the initial values
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
///////////////////////////////////////////////////////////////////////////////

#ifndef _INITIALIZATION_H
#define _INITIALIZATION_H

#include "data_structure.h"
#include "parameter_reader.h"
#include "sci_reader.h"
#include "utility.h"
#include "solver.h"
#include "geometry.h"


///////////////////////////////////////////////////////////////////////////////
/// Initialize the parameters
///
///\param para Pointer to FFD parameters
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int initialize(PARA_DATA *para);

///////////////////////////////////////////////////////////////////////////////
/// Set the default value for parameters
///
///\param para Pointer to FFD parameters
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void set_default_parameter(PARA_DATA *para);

///////////////////////////////////////////////////////////////////////////////
/// Set default initial values for simulation variables
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int set_initial_data (PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Loop through all the BC cells and find if tiles there
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int if_exist_tiles(PARA_DATA *para, REAL **var, int **BINDEX) ;

int init_para_simp(PARA_DATA *para, PARA_DATA_SIMP *para_simp);

//flatten var
int flat_var(PARA_DATA *para, REAL **var, REAL *var_flat);
// flatten bindex
int flat_index(PARA_DATA *para, int **BINDEX, int *bindex_flat);
int unflat_var(PARA_DATA *para, REAL **var, REAL *var_flat);

#endif
