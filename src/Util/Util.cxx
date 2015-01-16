#include "Util.h"

ClassImp(Util)

Util::Util(){
} 

const char* Util::GetRootFileNameFromPath (const char* path) const {
  string pathstr(path);
  int start = pathstr.find_last_of("/")+1; // cut directory path
  int end   = pathstr.length();
  const string fname = pathstr.substr(start,end);
  return fname.c_str();
}

