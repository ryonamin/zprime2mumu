#ifndef __UTIL_H__
#define __UTIL_H__

#include <string>
#include "TObject.h"
#include "TString.h"

using namespace std;

class Util : public TObject 
{
  public :
  Util();
  const char*   GetRootFileNameFromPath( const char* path ) const; 
  
  ClassDef(Util,1)
};
#endif
