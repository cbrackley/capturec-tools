//***************************************************************************
//
// Header for 
// Functions for working with bed files
//
//***************************************************************************

#ifndef BEDFILES_H
#define BEDFILES_H

#include<string>

using namespace std;

struct bedline {
  // data structure for bed file entry
  string chrom,
    name,
    strand;
  long int start,
    end;
  int istrand;
  bedline(const string&);
bedline() : chrom("0"), start(0), end(1) {};
  bool operator<(const bedline&) const;
};


struct bgdline {
  // data structure for bedGraph file entry
  string chrom;
  long int start,
    end;
  double value;
  bgdline(const string&);
};


#endif
