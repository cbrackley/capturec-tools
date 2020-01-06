//***************************************************************************
//
// Program to meaure the directionality from CaptureC data
//
//***************************************************************************

#include<iostream>
#include<cstdlib>
#include<string>
#include<map>
#include<fstream>
#include<sstream>
#include<cmath>

#include "directionality.h"
#include "bedfiles.h"

using namespace std;

int main(int argc, char *argv[]) {

  // get options from command line
  if (argc<7) {
    cout<<"Usage :"<<endl;
    cout<<"       ./directionality -t targetsfile -f inputslist -o outputfile -min MIN -max MAX"<<endl;
    cout<<"where       targetsfile  is a bed file with a list of targets."<<endl;
    cout<<"            inputslist   is a text file containing a list of paths to the normalized pile-up files along with the name of the target."<<endl;
    cout<<"            outfile      is a file name for the output."<<endl;
    cout<<"            MIN          is the minimum distance in bp from the target considered (Default 3000)."<<endl;
    cout<<"            MAX          is the maximum distance in bp from the target considered (Default 500,000)."<<endl;
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

  int min_dist=3000,
    max_dist=500000;

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
      min_dist = atoi(argv[argi+1]);
      argi += 2;

    } else if ( string(argv[argi]) == "-max" ) {
      // input file list
      if (!(argi+1 < argc)) {
        cerr<<"Error parsing command line (-max)"<<endl;
        exit(EXIT_FAILURE);
      }
      max_dist = atoi(argv[argi+1]);
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
  double dir;

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
  ouf<<"# chrom, start, end, targetname, directionality"<<endl;

  // Write messages
  cout<<"Finding directionalities for "<<inputfiles.size()<<" targets."<<endl;
  cout<<"Reads between "<<min_dist<<" and "<<max_dist<<" from the target centre are considered."<<endl;


  // Find directionalities
  for (map<string,string>::iterator it=inputfiles.begin();it!=inputfiles.end(); ++it) {
    inf.open( (it->second).c_str() );
    bool testfile=inf.good();
    inf.close();
    if ( !testfile ) {
      cerr<<" Warning : Cannot open file "<<it->second<<" skipping this."<<endl;
    } else {
      dir = get_directoinality(it->second,targets[it->first],max_dist,min_dist);
      ouf<<targets[it->first].chrom<<"\t"
	 <<targets[it->first].start<<"\t"
	 <<targets[it->first].end<<"\t"
	 <<it->first<<"\t"
	 <<dir<<endl;
    }

  }

  ouf.close();


}


double get_directoinality(const string &file,const bedline &trg,const int &max_dist,const int &min_dist) {
  // function to calculate the directionality
  ifstream inf;
  string line;
  double dir,
    upstream=0.0,
    downstream=0.0;
  int upwidth=0,
    downwidth=0;
  map<double,double> data;

  double trgmid=0.5*(trg.start+trg.end);
  double maxup=0,
    minup=10e9,
    maxdown=0,
    mindown=10e9;

  inf.open( file.c_str() );
  while ( getline(inf,line) ) {
    bgdline datapoint(line);

    if ( trg.chrom == datapoint.chrom &&
	 datapoint.end-trgmid>min_dist &&
	 datapoint.start-trgmid<max_dist ) {
	// downstream
	downstream += datapoint.value;
	if (datapoint.start<mindown) {mindown=datapoint.start;}
	if (datapoint.end>maxdown) {maxdown=datapoint.end;}
    }
    if ( trg.chrom == datapoint.chrom &&
	 datapoint.end-trgmid>-max_dist &&
	 datapoint.start-trgmid<-min_dist ) {
	// upstream
	upstream += datapoint.value;
	if (datapoint.start<minup) {minup=datapoint.start;}
	if (datapoint.end>maxup) {maxup=datapoint.end;}
    }

  }
  inf.close();

  // now find the log_2 ratio of the up/down stream reads per bp
  downwidth = maxdown-mindown;
  upwidth = maxup-minup;
  upstream /= double(upwidth);
  downstream /= double(downwidth);
  dir = log(upstream) - log(downstream);

  return dir;
}
