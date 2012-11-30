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

/**
 * standard constructor:
 */
TPSPnsM::TPSPnsM() : RecMat(TPTPnsM::gn(),PHM::gn()) { }

/**
 * copy constructor: constructs RecMat object
 * @param tpspm_c object that will be copied into this.
 */
TPSPnsM::TPSPnsM(const TPSPnsM &tpspm_c) : RecMat(tpspm_c){ }

/**
 * destructor
 */
TPSPnsM::~TPSPnsM(){ }

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:\n\n
 * @param a first sp index that forms the tp row index i together with b, part of the TPSPnsM row index
 * @param b second sp index that forms the tp row index i together with a, part of the TPSPnsM row index
 * @param c first sp index that forms the tp column index j together with d, part of the TPSPnsM row index
 * @param d second sp index that forms the tp column index j together with c, part of the TPSPnsM row index
 * @param e first sp index that forms part of the TPSPnsM colum index together with b
 * @param z second sp index that forms part of the TPSPnsM colum index together with a
 * @return the number on place TPSPnsM(i,j) with the right phase.
 */
double TPSPnsM::operator()(int a,int b,int c,int d,int e,int z) const{

   if( (a == b) || (c == d) )
      return 0;
   else{

      int phase = 1;

      if(a > b)
         phase *= -1;

      if(c > d)
         phase *= -1;

      int I = TPM::gs2t(a,b);
      int J = TPM::gs2t(c,d);

      int i = TPTPnsM::gt2tpmm(I,J);

      int j = PHM::gs2ph(e,z);

      return phase * (*this)(i,j);

   }

}

/**
 * construct a TPSPnsM object by reordering the indices of a DPM object, needed for the efficient
 * evaluation of the dirprodtrace() function.
 * @param dpm input DPM
 */
void TPSPnsM::reorder(const DPM &dpm){

   int I,J;

   int a,b,c,d,e,z;

   for(int i = 0;i < gm();++i){

      I = TPTPnsM::gtpmm2t(i,0);
      J = TPTPnsM::gtpmm2t(i,1);

      a = TPM::gt2s(I,0);
      b = TPM::gt2s(I,1);

      c = TPM::gt2s(J,0);
      d = TPM::gt2s(J,1);

      for(int j = 0;j < gn();++j){

         e = PHM::gph2s(j,0);
         z = PHM::gph2s(j,1);

         (*this)(i,j) = dpm(a,b,e,c,d,z);

      }
   }

}
