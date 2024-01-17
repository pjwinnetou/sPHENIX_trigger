#include <iostream>
#include "triggerEff.C"
#include "TROOT.h"
#include <vector>
#include "TEfficiency.h"
#include "set"

using namespace std;
int main(int argc, char *argv[])
{
  int argseed = atoi(argv[1]);
  TTree* tree=0;
  triggerEff t(tree, argseed);
  t.Loop();
  return 0;
}
