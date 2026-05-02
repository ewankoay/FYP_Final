///////////////////////////////////////////////////////////////////////////////
///
/// \file   main.c
///
/// \brief  main entrance of FFD program
///
/// \author Wei Tian
///         University of Miami, Schneider Electric
///         w.tian@umiami.edu, Wei.Tian@Schneider-Electric.com
///
/// \date   6/17/2017
///
///////////////////////////////////////////////////////////////////////////////

#include "ffd.h"

int main () {
  ffd_log ("main(): start FFD program with visualization capability", FFD_NEW);
  ffd(0);
}//end of main
