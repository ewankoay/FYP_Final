///////////////////////////////////////////////////////////////////////////////
///
/// \file   solver.c
///
/// \brief  Solver of FFD
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
#ifndef _SOLVER_H
#define _SOLVER_H

#include "data_structure.h"
#include "data_writer.h"
#include "diffusion.h"
#include "projection.h"
#include "advection.h"
#include "timing.h"
#include "solver_gs.h"
#include "boundary.h"
#include "utility.h"


///////////////////////////////////////////////////////////////////////////////
/// FFD solver
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int FFD_solver(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Calculate the temperature
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int temp_step(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Calculate the contaminant concentration
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int den_step(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Calculate the velocity
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int vel_step(PARA_DATA *para, REAL **var,int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Solver for equations
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param var_type Variable type
///\param Pointer to variable
///
///\return 0 if not error occurred
///////////////////////////////////////////////////////////////////////////////
int equ_solver(PARA_DATA *para, REAL **var, FFD_TERM which_term, int Type, REAL *x);

///////////////////////////////////////////////////////////////////////////////
/// Restore the conservation of scalar virables after advection using Semi-Lagrangian
/// method
/// Literature: https://engineering.purdue.edu/~yanchen/paper/2015-1.pdf
/// Wei Tian 6/20/2017, @ Schneider Electric, Andover, MA
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///\Param psi Pointer to scalar variables after advection
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int scalar_conservation(PARA_DATA *para, REAL **var, REAL *psi0, REAL *psi, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////////////
/// Check energy inconservation after advection
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///\Param psi Pointer to scalar variables after advection
///\Param psi0 Pointer to scalar variables before advection
///\return 0 if no error occurred
///\ Wei Tian
///\ 7/7/2017
//////////////////////////////////////////////////////////////////////////////////////
REAL adv_inconservation(PARA_DATA *para, REAL **var, REAL *psi0, REAL *psi, int **BINDEX);


int assign_tile_velocity(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////////////
/// Calculate the flow rates through the tiles using corrected pressure
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///\p_corr corrected pressure
///\return 0 if no error occurred
///\ Wei Tian
///\ 09/05/2017
//////////////////////////////////////////////////////////////////////////////////////
REAL pressure_correction(PARA_DATA *para, REAL **var, int **BINDEX, REAL p_corr);

///////////////////////////////////////////////////////////////////////////////
/// Check Imbalance (the unit is W, so there is no need to multiply a time step to inlet and outlet)
///
/// Wei Tian, update 6/20/2017, @Schneider Electric, Andover, MA
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int CheckImbalance(PARA_DATA *para, REAL **var, int var_type, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////////////
/// Assign the the pressure for the tiles after determining the mass flow rate proposed by Wei Tian
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///\return 0 if no error occurred
///\ Wei Tian
///\ 09/11/2017
//////////////////////////////////////////////////////////////////////////////////////
int assign_tile_velocity_hybrid(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////////////
/// dynamically update the pressure for the tiles after determining the mass flow rate proposed by Wei Tian
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///\return 0 if no error occurred
///\ Wei Tian
///\ 09/11/2017
//////////////////////////////////////////////////////////////////////////////////////
int update_tile_pressure(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////////////
///\Function Brief
/// Model the plenum and room as a whole. However, this is one of the two methods proposed in the issue 20.
/// The general idea of this method is that the tile is going to be treated as inlet cell for the above cell (in room),
/// and outlet cell for the below cell (in plenum). That being said, for the cells below the tile,
/// a Neumann BC is applid for velocity and a Dirichlet BC is applied for pressure. For the cells above
/// the tile, a Neumann BC is applied for the pressure and a Dirichlet BC is applied for the velocity.
/// Moreover, the pressure for one cell is always fixed at zero, so that the pressure is in a reasonable range.
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///\return 0 if no error occurred
///
///\ Wei Tian, 9-29-2017, Wei.Tian@Schneider-Electric.com
//////////////////////////////////////////////////////////////////////////////////////
int tile_room_split(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////////////
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///\return 0 if no error occurred
///
///\ Wei Tian, 10-19-2017, Wei.Tian@Schneider-Electric.com
//////////////////////////////////////////////////////////////////////////////////////
int tile_room_coupled(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////////////
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///\return 0 if no error occurred
///
///\ Wei Tian, 10-19-2017, Wei.Tian@Schneider-Electric.com
//////////////////////////////////////////////////////////////////////////////////////
int tile_source(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////////////
/// The black box model of rack, which treat the rack as a box with inlet outlet and heat dissipation
/// The temperature stratification of inlet temperature is kept in the outlet temperature
/// The velocity at inlet and outlet is the same
/// The inlet of rack is treated as outlet for the DC room while the outlet of rack is treated as inlet for DC room
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///\return 0 if no error occurred
///
///\ Wei Tian, 1-20-2018, Wei.Tian@Schneider-Electric.com
//////////////////////////////////////////////////////////////////////////////////////
int rack_model_black_box(PARA_DATA *para, REAL **var, int **BINDEX);

#endif
