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

#include "advection.h"
#include "boundary.h"
#include "data_structure.h"
#include "data_writer.h"
#include "diffusion.h"
#include "projection.h"
#include "solver_gs.h"
#include "timing.h"
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
int vel_step(PARA_DATA *para, REAL **var, int **BINDEX);

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
int equ_solver(PARA_DATA *para, REAL **var, FFD_TERM which_term, int Type,
               REAL *x);

///////////////////////////////////////////////////////////////////////////////
/// Restore the conservation of scalar virables after advection using
/// Semi-Lagrangian method Literature:
/// https://engineering.purdue.edu/~yanchen/paper/2015-1.pdf Wei Tian 6/20/2017,
/// @ Schneider Electric, Andover, MA
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///\Param psi Pointer to scalar variables after advection
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int scalar_conservation(PARA_DATA *para, REAL **var, REAL *psi0, REAL *psi,
                        int **BINDEX);

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
REAL adv_inconservation(PARA_DATA *para, REAL **var, REAL *psi0, REAL *psi,
                        int **BINDEX);

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
REAL pressure_correction(PARA_DATA *para, REAL **var, int **BINDEX,
                         REAL p_corr);

///////////////////////////////////////////////////////////////////////////////
/// Check Imbalance (the unit is W, so there is no need to multiply a time step
/// to inlet and outlet)
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
// calculate the initial tile flow rate based on flow network model
// Wei Tian
// 03/12/2019
//////////////////////////////////////////////////////////////////////////////////////
int initial_tile_velocity(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////////////
/// Assign the the velocity for the tiles after determining the pressure
/// correction A bisec method is used to solve the non-linear equations. For
/// more insights of the solver, refer to
/// http://cims.nyu.edu/~donev/Teaching/NMI-Fall2010/Lecture6.handout.pdf
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///\return 0 if no error occurred
///\ Wei Tian
///\ 09/05/2017
//////////////////////////////////////////////////////////////////////////////////////
int tile_pressure_correction_method(PARA_DATA *para, REAL **var, int **BINDEX);

int update_tile_pressure(PARA_DATA *para, REAL **var, int **BINDEX);

int update_tile_velocity(PARA_DATA *para, REAL **var, int **BINDEX);

int tile_room_split(PARA_DATA *para, REAL **var, int **BINDEX);

REAL flowrate_pressure_correction(PARA_DATA *para, REAL **var, int **BINDEX,
                                  REAL p_prime);
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
/// The black box model of rack, which treat the rack as a box with inlet outlet
/// and heat dissipation The temperature stratification of inlet temperature is
/// kept in the outlet temperature The velocity at inlet and outlet is the same
/// The inlet of rack is treated as outlet for the DC room while the outlet of
/// rack is treated as inlet for DC room
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///\return 0 if no error occurred
///
///\ Wei Tian, 1-20-2018, Wei.Tian@Schneider-Electric.com
//////////////////////////////////////////////////////////////////////////////////////
int rack_model_black_box(PARA_DATA *para, REAL **var, int **BINDEX);

#endif
