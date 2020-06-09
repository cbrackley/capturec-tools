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
  double midpoint();

private:
  double the_midpoint;
  bool set_midpoint;
};


struct bgdline {
  // data structure for bedGraph file entry
  string chrom;
  long int start,
    end;
  double value;
  bgdline(const string&);
  bgdline(const string &, const long int &, const long int &,const double &);
  bool operator<(const bgdline&) const;  
  double midpoint() const;

private:
  mutable double the_midpoint;
  mutable bool set_midpoint;

};




#endif
