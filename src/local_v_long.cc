//***************************************************************************
//
// Program to get the ratio between the number of local and long range 
// interactions, one value per target.
//
//***************************************************************************

#include<iostream>
#include<cstdlib>
#include<string>
#include<map>
#include<fstream>
#include<sstream>
#include<cmath>

#include "bedfiles.h"

using namespace std;

int main(int argc, char *argv[]) {

  // get options from command line
  if (argc<7) {
    cout<<"Usage :"<<endl;
    cout<<"       ./directionality -t targetsfile -f inputslist -o outputfile [-min MIN] [-max MAX] [-h THRESH]"<<endl;
    cout<<"where       targetsfile  is a bed file with a list of targets."<<endl;
    cout<<"            inputslist   is a text file containing a list of paths to the normalized pile-up files along with the name of the target."<<endl;
    cout<<"            outfile      is a file name for the output."<<endl;
    cout<<"            MIN          OPTIONAL: interactions closer than this (bp) are ignored (Default=1000)"<<endl;
    cout<<"            MAX          OPTIONAL: interactions further than this (bp) are ignored (Default=10,000,000)"<<endl;
    cout<<"            THRESH       OPTIONAL: cut off for where local ends and long range starts (bp) (Default=100,000)"<<endl;
    cout<<endl;
    cout<<"The inputslist file must have the following format:"<<endl;
    cout<<"         /path/to/captured_normalizedpileup_probe1.bdg        nameoftarget1"<<endl;
    cout<<"         /path/to/captured_normalizedpileup_probe2.bdg        nameoftarget2"<<endl;
    cout<<"where the named target must be present in the targetsfile"<<endl;
    exit(EXIT_FAILURE);
  }

  string targetsfile,
    inputslist,
    outputfile;

  double min_dist=1000,
    max_dist=10000000,
    cutoff=100000;

  int argi=1;
  while (argi < argc) {

    if ( string(argv[argi]) == "-t" ) {
      // targets file
      if (!(argi+1 < argc)) {
        cerr<<"Error parsing command line (-t)"<<endl;
        exit(EXIT_FAILURE);
      }
      targetsfile = string(argv[argi+1]);
      argi += 2;

    } else if ( string(argv[argi]) == "-o" ) {
      // output file
      if (!(argi+1 < argc)) {
        cerr<<"Error parsing command line (-o)"<<endl;
        exit(EXIT_FAILURE);
      }
      outputfile = string(argv[argi+1]);
      argi += 2;

    } else if ( string(argv[argi]) == "-f" ) {
      // input file list
      if (!(argi+1 < argc)) {
        cerr<<"Error parsing command line (-f)"<<endl;
        exit(EXIT_FAILURE);
      }
      inputslist = string(argv[argi+1]);
      argi += 2;

    } else if ( string(argv[argi]) == "-min" ) {
      // input file list
      if (!(argi+1 < argc)) {
        cerr<<"Error parsing command line (-min)"<<endl;
        exit(EXIT_FAILURE);
      }
      min_dist = atof(argv[argi+1]);
      argi += 2;

    } else if ( string(argv[argi]) == "-max" ) {
      // input file list
      if (!(argi+1 < argc)) {
        cerr<<"Error parsing command line (-max)"<<endl;
        exit(EXIT_FAILURE);
      }
      max_dist = atof(argv[argi+1]);
      argi += 2;

    } else if ( string(argv[argi]) == "-h" ) {
      // input file list
      if (!(argi+1 < argc)) {
        cerr<<"Error parsing command line (-h)"<<endl;
        exit(EXIT_FAILURE);
      }
      cutoff = atof(argv[argi+1]);
      argi += 2;

    } else {
      cerr<<"Error parsing command line (unrecognised option)"<<endl;
      exit(EXIT_FAILURE);
    }

  }

  // Set up variables
  ifstream inf;
  ofstream ouf;
  map<string,bedline> targets;
  map<string,string> inputfiles;

  string line;
  double sep,
    localCount=0.0,
    longCount=0.0,
    locscale,
    lonscale;

  // Read the targets file
  inf.open( targetsfile.c_str() );
  if ( !inf.good() ) {
    cerr<<" ERROR : Cannot open file "<<targetsfile<<endl;
    exit(EXIT_FAILURE);
  }
  while ( getline(inf,line) ) {
    bedline newtrg(line);
    targets[newtrg.name]=newtrg;
  }
  inf.close();


  // Read the inputs list
  inf.open( inputslist.c_str() );
  if ( !inf.good() ) {
    cerr<<" ERROR : Cannot open file "<<inputslist<<endl;
    exit(EXIT_FAILURE);
  }
  while ( getline(inf,line) ) {
    istringstream sline(line);
    string filename,
      trgname;
    sline>>filename>>trgname;
    inputfiles[trgname]=filename;
    if ( targets.count(trgname)==0 ) {
      cerr<<" ERROR : target "<<trgname<<" is not in the targets file."<<endl;
      exit(EXIT_FAILURE);
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
  ouf<<"# chrom, start, end, targetname, loc/long ratio, total local, total long"<<endl;
  ouf<<"# Interactions with regions between "<<min_dist<<" and "<<cutoff<<" bp are local"<<endl;
  ouf<<"# Interactions with regions between "<<cutoff<<" and "<<max_dist<<" bp are long range"<<endl;

  // Write messages
  cout<<"Finding local and long range reads for "<<inputfiles.size()<<" targets."<<endl;
  cout<<"Interactions with regions between "<<min_dist<<" and "<<cutoff<<" bp are local"<<endl;
  cout<<"Interactions with regions between "<<cutoff<<" and "<<max_dist<<" bp are long range"<<endl;

  // Parse input files
  for (map<string,string>::iterator it=inputfiles.begin();it!=inputfiles.end(); ++it) {
    inf.open( (it->second).c_str() );
    if ( !inf.good() ) {
      cerr<<" ERROR : Cannot open file "<<it->second<<" execution terminated."<<endl;
      exit(EXIT_FAILURE);
    } 
    double actual_loc_low=1e12,
      actual_loc_hi=0,
      actual_lon_low=1e12,
      actual_lon_hi=0;
    while ( getline(inf,line) ) {
      bgdline datapoint(line);
      if ( targets[it->first].chrom == datapoint.chrom ) { 
	// only consider same chrom as target
	sep =  abs(datapoint.midpoint()-targets[it->first].midpoint());
	if (sep>=min_dist && sep<=max_dist) {
	  if (sep<cutoff) {
	    // its local
	    localCount += datapoint.value;
	    if ( datapoint.midpoint()>targets[it->first].midpoint() ) {
	      if ( datapoint.start-targets[it->first].midpoint()<actual_loc_low ) {
		actual_loc_low=datapoint.start-targets[it->first].midpoint();
	      }
	      if ( datapoint.end-targets[it->first].midpoint()>actual_loc_hi ) {
		actual_loc_hi=datapoint.end-targets[it->first].midpoint();
	      }
	    } else {
	      if ( targets[it->first].midpoint()-datapoint.start>actual_loc_hi ) {
		actual_loc_hi=targets[it->first].midpoint()-datapoint.start;
	      }
	      if ( targets[it->first].midpoint()-datapoint.end>actual_loc_low ) {
		actual_loc_low=targets[it->first].midpoint()-datapoint.end;
	      }
	    }
	  } else {
	    // its long range
	    longCount += datapoint.value;
	    if ( datapoint.midpoint()>targets[it->first].midpoint() ) {
	      if ( datapoint.start-targets[it->first].midpoint()<actual_lon_low ) {
		actual_lon_low=datapoint.start-targets[it->first].midpoint();
	      }
	      if ( datapoint.end-targets[it->first].midpoint()>actual_lon_hi ) {
		actual_lon_hi=datapoint.end-targets[it->first].midpoint();
	      }
	    } else {
	      if ( targets[it->first].midpoint()-datapoint.start>actual_lon_hi ) {
		actual_lon_hi=targets[it->first].midpoint()-datapoint.start;
	      }
	      if ( targets[it->first].midpoint()-datapoint.end>actual_lon_low ) {
		actual_lon_low=targets[it->first].midpoint()-datapoint.end;
	      }
	    }
	  }
	}
      }
    }
    // adjustment factor to take into account different sized regions
    locscale=(actual_loc_hi-actual_loc_low)/(cutoff-min_dist);
    lonscale=(actual_lon_hi-actual_lon_low)/(max_dist-cutoff);
    //cout<<locscale<<" "<<lonscale<<endl;

    // output
    ouf<<targets[it->first].chrom<<"\t"
       <<targets[it->first].start<<"\t"
       <<targets[it->first].end<<"\t"
       <<it->first<<"\t"
       <<(localCount*locscale)/(longCount*lonscale)<<"\t"
       <<localCount*locscale<<"\t"
       <<longCount*lonscale<<endl;

    inf.close();
  }

  ouf.close();

}
