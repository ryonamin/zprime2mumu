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
  void setNSig( float in ) { fnsig = in; }
  private:
  TreeHandler& fT;
  float fnsig;
  ClassDef(MCParticleFinder,1)
};

