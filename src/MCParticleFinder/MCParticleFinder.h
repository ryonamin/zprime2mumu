#include "TObject.h"
#include "TTree.h"
#include <vector>

class TreeHandler;

class MCParticleFinder : public TObject
{
  public:
  MCParticleFinder(TreeHandler& in);
  ~MCParticleFinder();
  int getMatchedMCId(int recoId);
  void setNSig_Pt( float in )   { fnsig_pt = in; }
  void setNSig_Phi( float in )  { fnsig_phi = in; }
  void setNSig_Eta( float in )  { fnsig_eta = in; }
  void useChargeInfo( bool in ) { fuseChargeInfo = in; }
  private:
  TreeHandler& fT;
  float fnsig_pt;
  float fnsig_phi;
  float fnsig_eta;
  bool  fuseChargeInfo;
  ClassDef(MCParticleFinder,1)
};

