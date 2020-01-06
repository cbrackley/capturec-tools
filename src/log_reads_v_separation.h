//***************************************************************************
//
// Header for:
// Program to get mean reads vs genomic separation with logarythmicly 
// spaced bins.
//
//***************************************************************************

#ifndef LOGRDVSEP_H
#define LOGRDVSEP_H

#include<string>

#include "bedfiles.h"

using namespace std;

double get_log_bin(const double &,const double &);

double log_bin_width(const double &,const double &);

#endif
