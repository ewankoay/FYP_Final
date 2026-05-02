///////////////////////////////////////////////////////////////////////////////
///
/// \file   chen_zero_equ_model.h
///
/// \brief  Computes turbulent viscosity using Chen's zero equ model
///
/// \author Mingang Jin, Qingyan Chen
///         Purdue University
///         Jin55@purdue.edu, YanChen@purdue.edu
///         Wangda Zuo
///         University of Miami
///         W.Zuo@miami.edu
///
/// \date   8/3/2013
///
/// This file provides function that computes the turbulent viscosity using
/// Chen's zero equation model
///
///////////////////////////////////////////////////////////////////////////////
#ifndef _CHEN_ZERO_EQU_MODEL_H_
#define _CHEN_ZERO_EQU_MODEL_H_

#include "data_structure.h"


///////////////////////////////////////////////////////////////////////////////
/// Computes turbulent viscosity using Chen's zero equation model
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param i I-index of the control volume
///\param j J-index of the control volume
///\param k K-index of the control volume
///
///\return Turbulent Kinematic viscosity
///////////////////////////////////////////////////////////////////////////////
REAL nu_t_chen_zero_equ(PARA_DATA *para, REAL **var, int i, int j, int k);

///////////////////////////////////////////////////////////////////////////////
/// Computes turbulent thermal diffusivity using Chen's zero equation model
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\Author: Wei Tian @ Schneider Electric, Andover, MA
///
///\Initial Implementation: 7/11/2017
///\return Turbulent Kinematic viscosity
///
///\ This is to reimplement the turbulent thermal diffusivity for FFD using Chen's
///  zero equation model. Theoretically, by diving the turbulent viscosity over the
///  Pr number, we can easily get the turbulent themal diffusivity coefficients. However,
///  according to our test, that could cause severe energy imblance for the whole space.
///  The root reason to explain that is not fully identified yet. However, we presumbaly beleive that
///  the dradmatically (over 100) spatial differnce of the thermal diffusivity coefficiet may lead to
///  energy imbalance, happening in the diffusion term. We also tried using constant turbulence model,
///  which essentially multiply 100 to the lamninar one. Applying that in the FFD, we did not see
///  the deterioration of energy balance, at any level.
///\ Inspire by this, we propose to uniformly disseminate the turbulent thermal diffusivity over the
///  whole fluid space. According to our test, this method can perfectly achieves energy balance for
///  the room space, while keeping the turbulent features in some extent of the thermal environment.
///
///\ Last update: 7/11/2017
///////////////////////////////////////////////////////////////////////////////
REAL alpha_t_chen_zero_equ(PARA_DATA *para, REAL **var);

#endif
