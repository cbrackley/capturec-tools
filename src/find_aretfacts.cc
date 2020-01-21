//***************************************************************************
//
// Program to identify aretfacts between replicates
//
//***************************************************************************

#include<iostream>
#include<cstdlib>
#include<string>
#include<set>
#include<map>
#include<vector>
#include<fstream>
#include<sstream>
#include<cmath>
#include <sys/stat.h>


#include "find_aretfacts.h"
#include "bedfiles.h"

using namespace std;

bool DirExist(const std::string &s)
{
  struct stat buffer;
  return (stat (s.c_str(), &buffer) == 0);
}



int main(int argc, char *argv[]) {

  // get options from command line
  if (argc<9) {
    cout<<"Usage :"<<endl;
    cout<<"       ./find_aretfacts -f indir1 outdir1 -f indir2 outdir2 [-f indir3 outdir3...] -t target -a factor"<<endl;
    cout<<"where       indir1       is a directory with rep1 input files captured_rawpileup_ and captured_normalizedpileup_ files"<<endl;
    cout<<"            outdir1      is a directory for output of rep1."<<endl;
    cout<<"            target       is the name of a target."<<endl;
    cout<<"            a            is the factor for how much bigger than the other replicates the signal has to be to be condiered an aretfact."<<endl;
    cout<<endl;
    cout<<"The indirs and outdirs must come in pairs, and there must be at least two pairs."<<endl;
    exit(EXIT_FAILURE);
  }

  set<string> indir;
  map<string, string> outdir;
  string target;
  double factor=10.0;

  int argi=1;
  while (argi < argc) {

    if ( string(argv[argi]) == "-f" ) {
      // targets file
      if (!(argi+1 < argc)) {
        cerr<<"Error parsing command line (-f)"<<endl;
        exit(EXIT_FAILURE);
      }
      indir.insert( string(argv[argi+1]) );
      outdir[ string(argv[argi+1]) ] = string(argv[argi+2]);
      argi += 3;

    } else if ( string(argv[argi]) == "-t" ) {
      // targets file
      if (!(argi+1 < argc)) {
        cerr<<"Error parsing command line (-t)"<<endl;
        exit(EXIT_FAILURE);
      }
      target = string(argv[argi+1]);
      argi += 2;

    } else if ( string(argv[argi]) == "-a" ) {
      // input file list
      if (!(argi+1 < argc)) {
        cerr<<"Error parsing command line (-a)"<<endl;
        exit(EXIT_FAILURE);
      }
      factor = atof(argv[argi+1]);
      argi += 2;

    } else {
      cerr<<"Error parsing command line (unrecognised option)"<<endl;
      exit(EXIT_FAILURE);
    }

  }

  // test variables
  if ( indir.size()<2 ) {
      cerr<<"Not enough input dirs"<<endl;
      exit(EXIT_FAILURE);
  }
  if ( outdir.size()<2 ) {
      cerr<<"Not enough output dirs"<<endl;
      exit(EXIT_FAILURE);
  }
  if ( indir.size()!=outdir.size() ) {
      cerr<<"Must have same number of input and output dirs."<<endl;
      exit(EXIT_FAILURE);
  }
  if ( target=="" ) {
      cerr<<"Invalid target."<<endl;
      exit(EXIT_FAILURE);
  }
  if ( factor<=0 ) {
      cerr<<"Invalid value for factor."<<endl;
      exit(EXIT_FAILURE);
  }

  // write a message
  cout<<"Looking for artefacts in target "<<target<<" in directories"<<endl;
  for (set<string>::const_iterator it=indir.begin(); it!=indir.end(); ++it) {
    cout<<"   "<<*it<<"    --->   "<<outdir[*it]<<endl;
  }
  cout<<"Where an artefact is a read value greater than "<<factor<<" times the mean of the other valies"<<endl;

  // set up rest of variables
  ifstream inf;
  ofstream ouf;
  string line;

  map<string, set<bgdline> > input_rep;
  map<double,  datapoint> data;

  int counter=0;

  double sum,
    meanoftherest;

  // test output dir exists and files do not
  for (set<string>::const_iterator it=indir.begin(); it!=indir.end(); ++it) {
    if ( !DirExist(outdir[*it]) ) {
      cerr<<"Cannot find directory "<<outdir[*it]<<endl;
      exit(EXIT_FAILURE);
    }
    inf.open( (outdir[*it]+"/captured_rawpileup_"+target+".bdg").c_str() );
    if ( inf.good() ) {
      cerr<<"File "<<outdir[*it]+"/captured_rawpileup_"+target+".bdg"<<" already exists, will not overwrite."<<endl;
      exit(EXIT_FAILURE);
    }
    inf.close();
  }

  // load inputs into memory
  for (set<string>::const_iterator it=indir.begin(); it!=indir.end(); ++it) {
    //cout<<"Loading "<<*it+"/captured_normalizedpileup_"+target+".bdg"<<endl;
    inf.open( (*it+"/captured_normalizedpileup_"+target+".bdg").c_str() );
    if ( !inf.good() ) {
      cerr<<"Cannot open file "<<*it+"/captured_normalizedpileup_"+target+".bdg"<<endl;
      exit(EXIT_FAILURE);
    }
    input_rep[*it]=set<bgdline>();
    while ( getline(inf,line) ) {
      input_rep[*it].insert( bgdline(line) );
    }
    inf.close();
  }

  // convert the data into a list
  for (set<string>::const_iterator REPit=indir.begin(); REPit!=indir.end(); ++REPit) {
    for (set<bgdline>::iterator LINEit=input_rep[*REPit].begin(); LINEit!=input_rep[*REPit].end(); ++LINEit) {
      double mid=0.5*double(LINEit->start+LINEit->end);
      data[ mid ].chrom=LINEit->chrom;
      data[ mid ].start=LINEit->start;
      data[ mid ].end=LINEit->end;
      data[ mid ].value[*REPit]=LINEit->value;
    }
  }

  // save memory by clearing the sets
  input_rep.clear();

  // find artefacts
  for (map<double,  datapoint>::iterator it=data.begin(); it!=data.end(); ++it) {
    // fill zeros and get sum
    sum = 0.0;
    for (set<string>::const_iterator REPit=indir.begin(); REPit!=indir.end(); ++REPit) {
      if ( it->second.value.count( *REPit )==0 ) {
	it->second.value[*REPit] = 0;
      }
      sum += it->second.value[*REPit];
    }
    // test for aretefacts
    for (set<string>::const_iterator REPit=indir.begin(); REPit!=indir.end(); ++REPit) {
      meanoftherest = (sum-it->second.value[*REPit])/indir.size();
      if ( it->second.value[*REPit] > factor*meanoftherest ) {
	// it is an artefact
	it->second.artef[*REPit] = true;
	//cout<<it->second.value[*REPit] <<" "<< factor*meanoftherest <<endl;
	counter++;
      } else {
	// it is not
	it->second.artef[*REPit] = false;
      }
    }
  }
 
  cout<<"Found "<<counter<<" artefacts"<<endl;

  // copy input rawpileups into output rawpileups  
  for (set<string>::const_iterator it=indir.begin(); it!=indir.end(); ++it) {
    inf.open( (*it+"/captured_rawpileup_"+target+".bdg").c_str() );
    if ( !inf.good() ) {
      cerr<<"Cannot open file "<<*it+"/captured_rawpileup_"+target+".bdg"<<endl;
      exit(EXIT_FAILURE);
    }
    ouf.open( (outdir[*it]+"/captured_rawpileup_"+target+".bdg").c_str() );

    
    while ( getline(inf,line) ) {
      bgdline mybgdline = bgdline(line);
      if ( data[ mybgdline.midpoint() ].artef[*it] ) {
	// it is an artefact
	ouf<<data[ mybgdline.midpoint() ].chrom<<"\t"
	   <<data[ mybgdline.midpoint() ].start<<"\t"
	   <<data[ mybgdline.midpoint() ].end<<"\t"
	   <<"0"<<endl;
      } else {
	// it is not an artefact
	ouf<<line<<endl;
      }

    }

    ouf.close();
    inf.close();
  }


}
