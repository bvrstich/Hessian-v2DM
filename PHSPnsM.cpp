#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

vector< vector<int> > PHSPnsM::phmm2ph;
int **PHSPnsM::ph2phmm;

/**
 * initialize the static lists
 */
void PHSPnsM::init(){

  ph2phmm = new int * [PHM::gn()];

   for(int i = 0;i < PHM::gn();++i)
      ph2phmm[i] = new int [PHM::gn()];

   vector<int> v(2);

   int phmm = 0;

   //no symmetry, difference with t2hess
   for(int i = 0;i < PHM::gn();++i)
      for(int j = 0;j < PHM::gn();++j){

         v[0] = i;
         v[1] = j;

         phmm2ph.push_back(v);

         ph2phmm[i][j] = phmm;

         ++phmm;

      }

}

/**
 * deallocate the static lists
 */
void PHSPnsM::clear(){

   for(int i = 0;i < TPM::gn();++i)
      delete [] ph2phmm[i];

   delete [] ph2phmm;

}

/**
 * standard constructor:
 */
PHSPnsM::PHSPnsM() : RecMat(phmm2ph.size(),PHM::gn()) { }

/**
 * copy constructor: constructs RecMat object
 * @param phspm_c object that will be copied into this.
 */
PHSPnsM::PHSPnsM(const PHSPnsM &phspm_c) : RecMat(phspm_c){ }

/**
 * destructor
 */
PHSPnsM::~PHSPnsM(){ }

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:\n\n
 * @param a first sp index that forms the ph row index i together with b, part of the PHSPnsM row index
 * @param b second sp index that forms the ph row index i together with a, part of the PHSPnsM row index
 * @param c first sp index that forms the ph column index j together with d, part of the PHSPnsM row index
 * @param d second sp index that forms the ph column index j together with c, part of the PHSPnsM row index
 * @param e first sp index that forms part of the PHSPnsM colum index together with b
 * @param z second sp index that forms part of the PHSPnsM colum index together with a
 * @return the number on place PHSPnsM(i,j) with the right phase.
 */
double PHSPnsM::operator()(int a,int b,int c,int d,int e,int z) const{

   int I = PHM::gs2ph(a,b);
   int J = PHM::gs2ph(c,d);

   int i = ph2phmm[I][J];

   int j = PHM::gs2ph(e,z);

   return (*this)(i,j);

}

/**
 * construct a PHSPnsM object by reordering the indices of a PPHM object, needed for the efficient
 * evaluation of the dirprodtrace() function.
 * @param pphm input PPHM
 */
void PHSPnsM::reorder(const PPHM &pphm){

   int I,J;

   int a,b,c,d,e,z;

   for(int i = 0;i < gm();++i){

      I = phmm2ph[i][0];
      J = phmm2ph[i][1];

      a = PHM::gph2s(I,0);
      b = PHM::gph2s(I,1);

      c = PHM::gph2s(J,0);
      d = PHM::gph2s(J,1);

      for(int j = 0;j < gn();++j){

         e = PHM::gph2s(j,0);
         z = PHM::gph2s(j,1);

         (*this)(i,j) = pphm(e,a,b,z,c,d);

      }
   }

}
