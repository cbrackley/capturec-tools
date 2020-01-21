//***************************************************************************
//
// Header for 
// Program to identify aretfacts between replicates
//
//***************************************************************************

#ifndef FINDARTE_H
#define FINDARTE_H

#include<map>
#include<string>

using namespace std;

struct datapoint {

  double mid;
  string chrom;
  int  start,
    end;

  map<string, double> value;
  map<string, bool> artef;

};

#endif
