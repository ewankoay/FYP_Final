///////////////////////////////////////////////////////////////////////////////
///
/// \file   ffd_data_reader.c
///
/// \brief  Read the previous FFD result file (Tecplot format)
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

#ifndef _FFD_DATA_READER_H
#define _FFD_DATA_READER_H

#include "data_structure.h"

#include "utility.h"

FILE *file_old_ffd;

///////////////////////////////////////////////////////////////////////////////
/// Read the previous FFD simulation data in a format of standard output
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int read_ffd_data(PARA_DATA *para, REAL **var);

#endif
