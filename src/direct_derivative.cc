//***************************************************************************
//
// Program to calculate the derivative of a directionality profile
//
// Uses simplest linear interpolation.
// Does not care how far appart probes are.
//
//***************************************************************************

#include<iostream>
#include<cstdlib>
#include<string>
#include<set>
#include<fstream>
#include<sstream>

#include "bedfiles.h"

using namespace std;

int main(int argc, char *argv[]) {

  // get options from command line
  if (argc<5) {
    cout<<"Usage :"<<endl;
    cout<<"       ./direct_derivative -d directionalityfile -o outputfile"<<endl;
    cout<<"where       directionalityfile  is a ."<<endl;
    cout<<"            outfile      is a file name for the output."<<endl;
    cout<<endl;
    exit(EXIT_FAILURE);
  }

  string dirfile,
    inputslist,
    outputfile;

  int argi=1;
  while (argi < argc) {

    if ( string(argv[argi]) == "-d" ) {
      // targets file
      if (!(argi+1 < argc)) {
        cerr<<"Error parsing command line (-d)"<<endl;
        exit(EXIT_FAILURE);
      }
      dirfile = string(argv[argi+1]);
      argi += 2;
    
    } else if ( string(argv[argi]) == "-o" ) {
      // output file
      if (!(argi+1 < argc)) {
        cerr<<"Error parsing command line (-o)"<<endl;
        exit(EXIT_FAILURE);
      }
      outputfile = string(argv[argi+1]);
      argi += 2;

    } else {
      cerr<<"Error parsing command line (unrecognised option)"<<endl;
      exit(EXIT_FAILURE);
    }

  }

  // Set up variables
  ifstream inf;
  ofstream ouf;
  set<bgdline> targets;
  typedef set<bgdline>::const_iterator targit;

  string line;
  stringstream sline;
    
  
  // Read the dir file
  inf.open( dirfile.c_str() );
  if ( !inf.good() ) {
    cerr<<" ERROR : Cannot open file "<<dirfile<<endl;
    exit(EXIT_FAILURE);
  }
  while ( getline(inf,line) ) {
    if ( line.compare(0,1,"#")!=0 ) {
      string chrom,junk;
      long int start, end;
      double dir;
      sline.clear(); sline.str(line);
      sline>>chrom>>start>>end>>junk>>dir;
      targets.insert( bgdline(chrom,start,end,dir)  );
    }
  }
  inf.close();
  
  // Test output file, then open it
  inf.open( outputfile.c_str() );
  if ( inf.good() ) {
    cerr<<" ERROR : File "<<outputfile<<" already exists. Will not overwrite."<<endl;
    exit(EXIT_FAILURE);
  }
  inf.close();
  ouf.open( outputfile.c_str() );
  ouf<<"# chrom, start, end, derivative"<<endl;

  
  
  // Parse the targets in pairs
  double dx,dD,newpos;
  targit it1=targets.begin(),
    it2=++targets.begin();
  for ( ; it2!=targets.end(); ++it1,++it2 ) {
    
    if (it1->chrom != it2->chrom) {
      cerr<<" ERROR : Not all directionality entries are on the same chromosome."<<endl;
      cerr<<it1->chrom<<" "<<it2->chrom<<endl;
      exit(EXIT_FAILURE);
    }
   
    dx = it2->midpoint() - it1->midpoint();
    dD = it2->value - it1->value;
    newpos = it1->start + 0.5*dx;
    ouf<<it1->chrom<<"\t"
	<<int(newpos)<<"\t"
	<<int(newpos)+1<<"\t"
	<<dD/dx<<endl;
  }

  ouf.close();
  
}
