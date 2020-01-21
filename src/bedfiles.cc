//***************************************************************************
//
// Functions for working with bed files
//
//***************************************************************************

#include<string>
#include<sstream>

#include "bedfiles.h"

using namespace std;


bedline::bedline(const string &line) : set_midpoint(false) {
  // Convert a string into a bedline object
  stringstream sline;
  sline.str(line);
  sline>>chrom>>start>>end>>name>>strand;
  if (strand == "+") {
    istrand=1;
  } else if (strand == "-") {
    istrand=-1;
  } else {
    istrand=0;
  }
}



bool bedline::operator<(const bedline &b) const {
  if ( chrom != b.chrom ) {
    return chrom<b.chrom;
  } else if (start != b.start) {
    return start<b.start;
  } else {
    return end<b.end;
  }
}

double bedline::midpoint() {
  // Find the midpoint (bp) of the region
  if (!set_midpoint) {
    the_midpoint = 0.5*(start+end);
    set_midpoint = true;
  }
  return the_midpoint;
}



bgdline::bgdline(const string &line) : set_midpoint(false) {
  // Convert a string into a bedline object
  stringstream sline;
  sline.str(line);
  sline>>chrom>>start>>end>>value;
}

bool bgdline::operator<(const bgdline &b) const {
  if ( chrom != b.chrom ) {
    return chrom<b.chrom;
  } else if (start != b.start) {
    return start<b.start;
  } else {
    return end<b.end;
  }
}

double bgdline::midpoint() {
  // Find the midpoint (bp) of the region
  if (!set_midpoint) {
    the_midpoint = 0.5*(start+end);
    set_midpoint = true;
  }
  return the_midpoint;
}
