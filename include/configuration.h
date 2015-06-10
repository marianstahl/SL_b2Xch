/*

    2015-05-26 M. Stahl

    Class to dump stuff into and get out of. As this thing is growing I should think about doing this properly some day...

*/

#ifndef CONFIGURATION
#define CONFIGURATION
#include "TString.h"
#include "Riostream.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSystem.h"
#include <sys/utsname.h>
#include <vector>
#include <exception>
using namespace std;
enum Verbosity : int{quiet, normal, verbose, dump_all};
class configuration{

public:

  configuration(){
    debug = false;
    verb_lvl = normal;
    version = "v1";
    particle = "Xic";
    cutstring = "p_ProbNNp#GT#0,0.2:10";
    struct utsname buffer;uname(&buffer);
    environment(static_cast<TString>(buffer.nodename));
    nbinsforfit = 79;
    XcMlo = 2388.0; XcMhi = 2546.0;
  }

  configuration(const TString jobname){

    debug = false;
    verb_lvl = normal;

    try {TObjArray *djn = jobname.Tokenize("/"); if(djn->GetEntries() != 3) throw 1;}
    catch (int e){cout << "Exception: " << e << "\n Using default configuration. Please specify version/particle/cutstring" << endl; configuration();}
    version = (TString)((TObjString*)(jobname.Tokenize("/")->At(0)))->String();
    particle = (TString)((TObjString*)(jobname.Tokenize("/")->At(1)))->String();
    cutstring = (TString)((TObjString*)(jobname.Tokenize("/")->At(2)))->String();

    struct utsname buffer;uname(&buffer);
    environment(static_cast<TString>(buffer.nodename));

    nbinsforfit = 79;
    if(particle.Contains("Omegac")){XcMlo = 2616.0; XcMhi = 2774.0;}
    else {XcMlo = 2388.0; XcMhi = 2546.0;}
  }

  ~configuration(){}

  //parse machine (so far only heidelberg and cern machines)
  void environment(TString machine){
    printf("You are working on %s\n",machine.Data());
    if(machine.Contains("eidelberg") || machine.Contains("sigma0")){
      lxplus = false;
      dumpdir = "/auto/data/mstahl/dump_area/SLBaryonSpectroscopy/"+version;
      notedir = "~mstahl/Documents/AnaNotes/SLBaryonSpectroscopy/latest/latex";
      tupledir = "/auto/data/mstahl/SLBaryonSpectroscopy/";
    }
    else{
      lxplus = true;
      dumpdir = "/afs/cern.ch/work/m/mstahl/public/dumpdir/SLBaryonSpectroscopy/"+version;
      notedir = "~mstahl/Documents/AnaNotes/SLBaryonSpectroscopy/latest/latex";
      tupledir = "/eos/lhcb/users/m/mstahl/SLBaryonSpectroscopy/";
    }
    if(!gSystem->OpenDirectory(dumpdir))gSystem->mkdir(dumpdir,true);
    if(!gSystem->OpenDirectory(notedir))gSystem->mkdir(notedir,true);
    if(!gSystem->OpenDirectory(tupledir)){cout << "tuple-directory not found. terminating..." << endl; terminate();}
  }
  //simple getters/setters  
  TString get_particle(){return particle;}
  TString get_version(){return version;}
  TString get_cutstring(){return cutstring;}
  bool is_debug(){return debug;}
  Verbosity get_verbosity(){return verb_lvl;}
  TString get_dumpdir(){return dumpdir;}
  TString get_notedir(){return notedir;}
  TString get_tupledir(){return tupledir;}

  void verbosity(const int lvl, const bool dbg = 0){
    verb_lvl = static_cast<Verbosity>(lvl);
    debug = dbg;return;
  }
  void verbosity(Verbosity lvl, const bool dbg = 0){
    verb_lvl = lvl;
    debug = dbg;return;
  }

  void set_particle(const TString part){
    if(!(part.Contains("Xic") || part.Contains("Omegac")))cout << "sure about the particle?" << endl;
      particle = part;
      if(particle.Contains("Omegac")){XcMlo = 2616.0; XcMhi = 2774.0;}
      else {XcMlo = 2388.0; XcMhi = 2546.0;}
      return;
  }
  void set_version(const int ver){version.Form("v%d",ver);return;}
  void set_version(const TString ver){version = ver;return;}

  void fill_cs(TString cs){cs_vec.push_back(cs);return;}
  vector<TString> get_cs_vec(){return cs_vec;}
  TString get_cs_at(unsigned int i){
    try{return cs_vec.at(i);}
    catch(exception e){
      cout << e.what() << ", terminating because trying to access cutstring in cs_vec that doesn't exist..." << endl;
      terminate();
    }
  }
  void set_current_cs(TString cs){current_cs = cs;return;}
  TString get_current_cs(){return current_cs;}

  double get_XcMlo(){return XcMlo;}
  double get_XcMhi(){return XcMhi;}
  int get_NBinsForFit(){return nbinsforfit;}

private:

  TString particle;
  TString version;
  TString cutstring;
  vector<TString> cs_vec;
  TString current_cs;  

  bool debug;
  Verbosity verb_lvl;
  bool lxplus;

  TString dumpdir;
  TString notedir;
  TString tupledir;

  int nbinsforfit;
  double XcMlo;
  double XcMhi;

};
#endif
