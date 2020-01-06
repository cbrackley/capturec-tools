//***************************************************************************
//
// Program to get mean reads vs genomic separation with logarythmicly 
// spaced bins.
//
//***************************************************************************

#include<iostream>
#include<cstdlib>
#include<string>
#include<map>
#include<fstream>
#include<sstream>
#include<cmath>

#include "log_reads_v_separation.h"
#include "bedfiles.h"

using namespace std;

int main(int argc, char *argv[]) {

  // get options from command line
  if (argc<7) {
    cout<<"Usage :"<<endl;
    cout<<"       ./directionality -t targetsfile -f inputslist -o outputfile [-b LBW]"<<endl;
    cout<<"where       targetsfile  is a bed file with a list of targets."<<endl;
    cout<<"            inputslist   is a text file containing a list of paths to the normalized pile-up files along with the name of the target."<<endl;
    cout<<"            outfile      is a file name for the output."<<endl;
    cout<<"            LBW          OPTIONAL: logarythmic bin width (Default=0.25)"<<endl;
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

  double log_binwidth=0.25;

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

    } else if ( string(argv[argi]) == "-b" ) {
      // input file list
      if (!(argi+1 < argc)) {
        cerr<<"Error parsing command line (-b)"<<endl;
        exit(EXIT_FAILURE);
      }
      log_binwidth = atof(argv[argi+1]);
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
  double bincentre,
    actual_binwidth;
  map<double,double> sumReads,
    sum2Reads,
    errorReads,
    counterReads;


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


  // Write messages
  cout<<"Finding reads vs genomic separation, averging over "<<inputfiles.size()<<" targets."<<endl;
  cout<<"Using logarythmic bin width of "<<log_binwidth<<endl;


  // Parse input files
  for (map<string,string>::iterator it=inputfiles.begin();it!=inputfiles.end(); ++it) {
    inf.open( (it->second).c_str() );
    if ( !inf.good() ) {
      cerr<<" ERROR : Cannot open file "<<it->second<<" execution terminated."<<endl;
      exit(EXIT_FAILURE);
    } 
    while ( getline(inf,line) ) {
      bgdline datapoint(line);
      if ( targets[it->first].chrom == datapoint.chrom ) { 
	// only consider same chrom as target
	// normalize by bin width
	bincentre = get_log_bin( abs(datapoint.midpoint()-targets[it->first].midpoint()) ,log_binwidth );
	actual_binwidth = log_bin_width(bincentre,log_binwidth);
	sumReads[bincentre] += datapoint.value/actual_binwidth;
	sum2Reads[bincentre] += datapoint.value*datapoint.value/actual_binwidth/actual_binwidth;
	counterReads[bincentre] ++;
      }
    }
    inf.close();
  }

  // Finish averages
  for (map<double,double>::iterator it=sumReads.begin(); it!=sumReads.end(); ++it) {
    sumReads[it->first] /= double(counterReads[it->first]);
    sum2Reads[it->first] /= double(counterReads[it->first]);
    errorReads[it->first] = sqrt(sum2Reads[it->first]-sumReads[it->first]*sumReads[it->first])/sqrt(counterReads[it->first]);
  }

  // output
  ouf.open( outputfile.c_str() );
  ouf<<"# Log[bin cenre], Log[reads per kbp], standard error (log), bin centre, reads per kbp"<<endl;
  for (map<double,double>::iterator it=sumReads.begin(); it!=sumReads.end(); ++it) {
    ouf<<log(it->first)<<" "
       <<log(it->second*1000)<<" "
       <<errorReads[it->first]/it->second<<" "
       <<it->first<<" "
       <<it->second*1000<<" "
       <<endl;
  }
  ouf.close();


}



double get_log_bin(const double &pos,const double &logwidth ) {
  // convert a position to the centre of a logarythmically spaced bin
  double logpos=log(pos);
  return int( exp(int(logpos/logwidth)*logwidth+0.5*logwidth) );
}

double log_bin_width(const double &centre,const double &logwidth ) {
  // From a bin centre, return the bin width in bp
  double lcentre=log(centre);
  return exp(lcentre+0.5*logwidth) - exp(lcentre-0.5*logwidth);
}
