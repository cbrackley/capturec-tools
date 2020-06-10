//***************************************************************************
//
// Program to calculate some statistics from a set of binned profiles.
//
// Looks at the first 30Mbp of the chromsome (0 to 30Mbp)
// Might be useful for comparisons between replicates.
//
// Reads binned profiles. The larger the bin size the better
// (all profiles should have the same bins)
//
//***************************************************************************

#include<iostream>
#include<cstdlib>
#include<cmath>
#include<string>
#include<set>
#include<map>
#include<vector>
#include<fstream>
#include<sstream>
#include <algorithm>

#include "bedfiles.h"

#define MAXREGION 30000000

using namespace std;

struct Lstats {

  double loWisk,
    Q1,
    median,
    Q3,
    hiWisk;

  Lstats (const double &lW,const double &qq1,
	  const double &me,const double &qq3,const double &hW) :
    loWisk(lW), Q1(qq1), median(me), Q3(qq3), hiWisk(hW) {}
  
};

pair<double,double> get_prpnAtoB(const string&,const bedline&,const int&,const int&);
Lstats get_Lest0to30M(const string&,const bedline&);

int main(int argc, char *argv[]) {

  // get options from command line
  if (argc<7) {
    cout<<"Usage :"<<endl;
    cout<<"       ./read_starts -t targetsfile -f inputslist -o outputfile"<<endl;
    cout<<"where       targetsfile  is a bed file with a list of targets."<<endl;
    cout<<"            inputslist   is a text file containing a list of paths to the normalized pile-up files along"<<endl
	<<"                         with the name of the target and the name of the condition/replicate."<<endl;
    cout<<"            outfile      is the first part of a file name for the output."<<endl;
    cout<<endl;
    cout<<"The inputslist file must have the following format:"<<endl;
    cout<<"         /path/to/captured_normalizedpileup_probe1.bdg  nameof_condition_or_replicate   nameoftarget1"<<endl;
    cout<<"         /path/to/captured_normalizedpileup_probe2.bdg  nameof_condition_or_replicate   nameoftarget2"<<endl;
    cout<<"where the named target must be present in the targetsfile"<<endl;
    exit(EXIT_FAILURE);
  }

  string targetsfile,
    inputslist,
    outputfilestart;

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
      outputfilestart = string(argv[argi+1]);
      argi += 2;

    } else if ( string(argv[argi]) == "-f" ) {
      // input file list
      if (!(argi+1 < argc)) {
        cerr<<"Error parsing command line (-f)"<<endl;
        exit(EXIT_FAILURE);
      }
      inputslist = string(argv[argi+1]);
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
  set<string> conditions;
  typedef set<string>::const_iterator cond_it;
  set<string> list_of_all_targets;
  typedef set<string>::const_iterator alltrg_it;
  map<string, map<string,string> > inputfiles;
  typedef map<string, map<string,string> >::const_iterator ifilesC_it;
  typedef map<string,string>::const_iterator ifilesT_it;

  string line;

  map<string, map<string,pair<double,double>> > prpn0to30M,
    prpn0to10M,
    prpn10to20M,
    prpn20to30M;
  map<string, map<string,Lstats> > stats;
  int ci;
  
  // Test output files do not exist
  inf.open( (outputfilestart+"prpn0to30M.dat").c_str() );
  if ( inf.good() ) {
    cerr<<" ERROR : File "<<outputfilestart+"prpn0to30M.dat"<<" already exists. Will not overwrite."<<endl;
    exit(EXIT_FAILURE);
  }
  inf.close();

  inf.open( (outputfilestart+"prpn0to10M.dat").c_str() );
  if ( inf.good() ) {
    cerr<<" ERROR : File "<<outputfilestart+"prpn0to10M.dat"<<" already exists. Will not overwrite."<<endl;
    exit(EXIT_FAILURE);
  }
  inf.close();

  inf.open( (outputfilestart+"prpn10to20M.dat").c_str() );
  if ( inf.good() ) {
    cerr<<" ERROR : File "<<outputfilestart+"prpn10to20M.dat"<<" already exists. Will not overwrite."<<endl;
    exit(EXIT_FAILURE);
  }
  inf.close();

  inf.open( (outputfilestart+"prpn20to30M.dat").c_str() );
  if ( inf.good() ) {
    cerr<<" ERROR : File "<<outputfilestart+"prpn20to30M.dat"<<" already exists. Will not overwrite."<<endl;
    exit(EXIT_FAILURE);
  }
  inf.close();

  
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
      condname,
      trgname;
    sline>>filename>>condname>>trgname;
    conditions.insert( condname );
    list_of_all_targets.insert( trgname );
    inputfiles[condname][trgname]=filename;
    if ( targets.count(trgname)==0 ) {
      cerr<<" ERROR : target "<<trgname<<" is not in the targets file."<<endl;
      exit(EXIT_FAILURE);
    }
    if ( targets[trgname].start < MAXREGION ) {
      cerr<<" WARNING : target "<<trgname<<" is within the first 30Mb of the chromosome. Results will not make sense."<<endl;
      //exit(EXIT_FAILURE);
    }
  }
  inf.close();


  // Write messages
  cout<<"Finding some statistics on "<<inputfiles.size()<<" conditions:"<<endl;
  for (ifilesC_it it=inputfiles.begin(); it!=inputfiles.end(); ++it) {
    cout<<"  "<<it->first<<" with "<<it->second.size()<<" targets"<<endl;
  }


  // Loop round conditions
  for (cond_it cond=conditions.begin(); cond!=conditions.end(); ++cond) {

    // Loop round targets in this condition
    for (ifilesT_it trg=inputfiles[*cond].begin(); trg!=inputfiles[*cond].end(); ++trg) {
      prpn0to30M[*cond][trg->first] = get_prpnAtoB(trg->second,targets[trg->first],0,30000000);
      prpn0to10M[*cond][trg->first] = get_prpnAtoB(trg->second,targets[trg->first],0,10000000);
      prpn10to20M[*cond][trg->first] = get_prpnAtoB(trg->second,targets[trg->first],10000000,20000000);
      prpn20to30M[*cond][trg->first] = get_prpnAtoB(trg->second,targets[trg->first],20000000,30000000);

      stats[*cond].insert( pair<string,Lstats>( trg->first,get_Lest0to30M(trg->second,targets[trg->first]) ) );
    }
    
  }

  // Ouput prpn0to30M
  ouf.open( (outputfilestart+"prpn0to30M.dat").c_str() );
  ouf<<"# propotion of reads (and error) within the first 30Mb for all targets and conditions"<<endl;
  ouf<<"# Column 1: name of target"<<endl;
  ci=2;
  for (cond_it cond=conditions.begin(); cond!=conditions.end(); ++cond) {
    ouf<<"# Columns "<<ci<<" and "<<ci+1<<": proportion and error for condition "<<*cond<<endl;
    ci+=2;
  }
  for (alltrg_it trg=list_of_all_targets.begin(); trg!=list_of_all_targets.end(); ++trg) {
    ouf<<*trg;
    for (cond_it cond=conditions.begin(); cond!=conditions.end(); ++cond) {
      if ( prpn0to30M[*cond].find(*trg)==prpn0to30M[*cond].end() ) {
	// The target was not present for this condition
	ouf<<" "<<"0 0";
      } else {
	ouf<<" "<<prpn0to30M[*cond][*trg].first<<" "<<prpn0to30M[*cond][*trg].second;
      }
    }
    ouf<<endl;
  }
  ouf.close();


  // Ouput prpn0to10M
  ouf.open( (outputfilestart+"prpn0to10M.dat").c_str() );
  ouf<<"# propotion of reads (and error) within the first 30Mb for all targets and conditions"<<endl;
  ouf<<"# Column 1: name of target"<<endl;
  ci=2;
  for (cond_it cond=conditions.begin(); cond!=conditions.end(); ++cond) {
    ouf<<"# Columns "<<ci<<" and "<<ci+1<<": proportion and error for condition "<<*cond<<endl;
    ci+=2;
  }
  for (alltrg_it trg=list_of_all_targets.begin(); trg!=list_of_all_targets.end(); ++trg) {
    ouf<<*trg;
    for (cond_it cond=conditions.begin(); cond!=conditions.end(); ++cond) {
      if ( prpn0to10M[*cond].find(*trg)==prpn0to10M[*cond].end() ) {
	// The target was not present for this condition
	ouf<<" "<<"0 0";
      } else {
	ouf<<" "<<prpn0to10M[*cond][*trg].first<<" "<<prpn0to10M[*cond][*trg].second;
      }
    }
    ouf<<endl;
  }
  ouf.close();

  // Ouput prpn10to20M
  ouf.open( (outputfilestart+"prpn10to20M.dat").c_str() );
  ouf<<"# propotion of reads (and error) within the first 30Mb for all targets and conditions"<<endl;
  ouf<<"# Column 1: name of target"<<endl;
  ci=2;
  for (cond_it cond=conditions.begin(); cond!=conditions.end(); ++cond) {
    ouf<<"# Columns "<<ci<<" and "<<ci+1<<": proportion and error for condition "<<*cond<<endl;
    ci+=2;
  }
  for (alltrg_it trg=list_of_all_targets.begin(); trg!=list_of_all_targets.end(); ++trg) {
    ouf<<*trg;
    for (cond_it cond=conditions.begin(); cond!=conditions.end(); ++cond) {
      if ( prpn10to20M[*cond].find(*trg)==prpn10to20M[*cond].end() ) {
	// The target was not present for this condition
	ouf<<" "<<"0 0";
      } else {
	ouf<<" "<<prpn10to20M[*cond][*trg].first<<" "<<prpn10to20M[*cond][*trg].second;
      }
    }
    ouf<<endl;
  }
  ouf.close();


  // Ouput prpn20to30M
  ouf.open( (outputfilestart+"prpn20to30M.dat").c_str() );
  ouf<<"# propotion of reads (and error) within the first 30Mb for all targets and conditions"<<endl;
  ouf<<"# Column 1: name of target"<<endl;
  ci=2;
  for (cond_it cond=conditions.begin(); cond!=conditions.end(); ++cond) {
    ouf<<"# Columns "<<ci<<" and "<<ci+1<<": proportion and error for condition "<<*cond<<endl;
    ci+=2;
  }
  for (alltrg_it trg=list_of_all_targets.begin(); trg!=list_of_all_targets.end(); ++trg) {
    ouf<<*trg;
    for (cond_it cond=conditions.begin(); cond!=conditions.end(); ++cond) {
      if ( prpn20to30M[*cond].find(*trg)==prpn20to30M[*cond].end() ) {
	// The target was not present for this condition
	ouf<<" "<<"0 0";
      } else {
	ouf<<" "<<prpn20to30M[*cond][*trg].first<<" "<<prpn20to30M[*cond][*trg].second;
      }
    }
    ouf<<endl;
  }
  ouf.close();


  // Ouput boxplot
  ouf.open( (outputfilestart+"boxplotTo30M.dat").c_str() );
  ouf<<"# Groups of five columns give values to draw box plots (loWhisker, Q1, median, Q3, hiWhiser)"<<endl;
  ouf<<"# Column 1: name of target"<<endl;
  ci=2;
  for (cond_it cond=conditions.begin(); cond!=conditions.end(); ++cond) {
    ouf<<"# Columns "<<ci<<" to "<<ci+4<<": boxplot values for condition "<<*cond<<endl;
    ci+=5;
  }
  for (alltrg_it trg=list_of_all_targets.begin(); trg!=list_of_all_targets.end(); ++trg) {
    ouf<<*trg;
    for (cond_it cond=conditions.begin(); cond!=conditions.end(); ++cond) {
      if ( stats[*cond].find(*trg)==stats[*cond].end() )  {
	// The target was not present for this condition
	ouf<<" "<<"0 0 0 0 0";
      } else {
	ouf<<" "<<stats[*cond].find( *trg )->second.loWisk
	   <<" "<<stats[*cond].find( *trg )->second.Q1
	   <<" "<<stats[*cond].find( *trg )->second.median
	   <<" "<<stats[*cond].find( *trg )->second.Q3
	   <<" "<<stats[*cond].find( *trg )->second.hiWisk;
      }
    }
    ouf<<endl;
  }
  ouf.close();

  
}


pair<double,double> get_prpnAtoB(const string &file,const bedline &trg,const int &from,const int &to) {

  ifstream inf;
  string line;
  
  double sumTotal=0.0,
    sumFirst30=0.0;
  double prpn,
    error;
   
  inf.open( file.c_str() );
  while ( getline(inf,line) ) {
    bgdline datapoint(line);
    if ( trg.chrom == datapoint.chrom ) {
      sumTotal += datapoint.value;
      if ( datapoint.midpoint()>from && datapoint.midpoint()<=to ) {
	sumFirst30 += datapoint.value;
      }
    }
  }
  inf.close();

  prpn = sumFirst30/sumTotal;
  error = prpn*(1-prpn)/sqrt(sumTotal);
  
  return make_pair(prpn,error);

}


Lstats get_Lest0to30M(const string &file,const bedline &trg) {

  ifstream inf;
  string line;
  
  vector<double> A;
  long int region=30000000,
    delta_x=MAXREGION,
    last=-1,
    extra_zeros,
    counter=0;
  double sum_total=0.0,
    median,
    Q1,Q3,
    IQR15,
    loWisk,
    hiWisk;
  int n;
  
  inf.open( file.c_str() );
  while ( getline(inf,line) ) {
    bgdline datapoint(line);
    if ( trg.chrom == datapoint.chrom &&
	 datapoint.midpoint()<=region ) {
      A.push_back( datapoint.value );
      counter++;
    }
    if (last!=-1 && datapoint.start-last<delta_x) { // find the bin size
      delta_x = datapoint.start-last;
    }
    last = datapoint.start;
    sum_total += datapoint.value;
  }
  inf.close();

  // Now add in the missing zeros
  extra_zeros = int(double(region)/double(delta_x))-counter;
  for (int i=0;i<extra_zeros;i++) {
    A.push_back( 0 );
  }

  // Normalize by the total number of reads
  for (int i=0;i<A.size();i++) {
    A[i]/=sum_total;
  }
  
  // sort the vector
  sort( A.begin(), A.end() );
  
  // find median
  if ( A.size() % 2 == 0 ) {
    // even
    median = 0.5*(A[ int(0.5*(A.size())) -1 ] + A[ int(0.5*(A.size())) ]);  // need '-1' because counting is from 0
  } else {
    // odd
    median = A[ int(0.5*(A.size()+1)) -1 ];  // need '-1' because counting is from 0
  }

  // find Q1 and Q3
  if ( A.size() % 2 == 0 ) { // even, so look at first A.size()/2 values
    n=A.size()/2;
  } else { // odd so look at first (A.size()+1)/2 values
    n=(A.size()+1)/2;
  }
    
  if ( n % 2 == 0 ) {
    // even
    Q1 = 0.5*( A[ n/2 -1 ] + A [ n/2 ] ); // need '-1' because counting is from 0
    Q3 = 0.5*( A[ n + n/2 -1 ] + A [ n + n/2 ] ); // need '-1' because counting is from 0
  } else {
    // odd
    Q1 = A[ (n+1)/2 -1 ];
    Q3 = A[ n + (n+1)/2 -1 ];
  }

  // For the whiskers, take first and last values inside 1.5*IQR about the median
  IQR15 = 1.5*(Q3-Q1);
  loWisk=A[0];
  for (int i=0;A[i]<Q1-IQR15 && i+1<A.size();++i) {
    loWisk=A[i+1];  // lowest value inside Q1-IQR15
  }
  hiWisk=A.back();
  for (int i=0;A[i]<Q3+IQR15;++i) {
    hiWisk=A[i]; // highest vlue inside Q3+IQR15
  }

  Lstats stats(loWisk,Q1,median,Q3,hiWisk);
  return stats;
  
}
