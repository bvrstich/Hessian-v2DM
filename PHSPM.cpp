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

int **PHSPM::ph2phmm;
vector< vector<int> > PHSPM::phmm2ph;

/**
 * initialize the static lists
 */
void PHSPM::init(){

   ph2phmm = new int * [PHM::gn()];

   for(int i = 0;i < PHM::gn();++i)
      ph2phmm[i] = new int [PHM::gn()];

   vector<int> v(2);

   int phmm = 0;

   for(int i = 0;i < PHM::gn();++i)
      for(int j = i;j < PHM::gn();++j){

         v[0] = i;
         v[1] = j;

         phmm2ph.push_back(v);

         ph2phmm[i][j] = phmm;
         ph2phmm[j][i] = phmm;

         ++phmm;

      }

}

/**
 * deallocate the static lists
 */
void PHSPM::clear(){

   for(int i = 0;i < PHM::gn();++i)
      delete [] ph2phmm[i];

   delete [] ph2phmm;

}

/**
 * standard constructor:
 */
PHSPM::PHSPM() : RecMat(phmm2ph.size(),TPSPM::gspmmdim()) { }

/**
 * copy constructor: constructs RecMat object 
 * @param phspm_c object that will be copied into this.
 */
PHSPM::PHSPM(const PHSPM &phspm_c) : RecMat(phspm_c){ }

/**
 * destructor
 */
PHSPM::~PHSPM(){ }

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:\n\n
 * @param a first sp index that forms the tp row index i together with b, part of the PHSPM row index
 * @param b second sp index that forms the tp row index i together with a, part of the PHSPM row index
 * @param c first sp index that forms the tp column index j together with d, part of the PHSPM row index
 * @param d second sp index that forms the tp column index j together with c, part of the PHSPM row index
 * @param e first sp index that forms part of the PHSPM colum index together with b
 * @param z second sp index that forms part of the PHSPM colum index together with a
 * @return the number on place PHSPM(i,j)
 */
double PHSPM::operator()(int a,int b,int c,int d,int e,int z) const{

   int I = PHM::gs2ph(a,b);
   int J = PHM::gs2ph(c,d);

   int i = ph2phmm[I][J];

   int j = TPSPM::gs2spmm(e,z);

   return (*this)(i,j);

}

/**
 * access to the phmm list from outside of the class
 */
int PHSPM::gph2phmm(int I,int J){

   return ph2phmm[I][J];

}

/**
 * access to the phmm list from outside of the class
 * @param option == 0 return I, == 1 return J
 */
int PHSPM::gphmm2ph(int i,int option){

   return phmm2ph[i][option];

}

/**
 * fill the PHSPM object by tracing out one pair of indices of the symmetrized direct product of a PHM object with itsself
 * @param scale factor to scale the trace with
 */
void PHSPM::dirprodtrace(double scale,const PHM &G){

   int I,J;

   int a,b,c,d,e,z;

   for(int i = 0;i < gm();++i){//loops over PHMM space

      I = phmm2ph[i][0];
      J = phmm2ph[i][1];

      a = PHM::gph2s(I,0);
      b = PHM::gph2s(I,1);

      c = PHM::gph2s(J,0);
      d = PHM::gph2s(J,1);

      for(int j = 0;j < gn();++j){//loops over SPMM space

         e = TPSPM::gspmm2s(j,0);
         z = TPSPM::gspmm2s(j,1);

         (*this)(i,j) = 0.0;

         for(int l = 0;l < Tools::gM();++l)
            (*this)(i,j) += G(a,b,e,l) * G(c,d,z,l) + G(c,d,e,l) * G(a,b,z,l);

         (*this)(i,j) *= scale;

      }
   }

}

/**
 * @return the dimension of particle-hole matrix space
 */
int PHSPM::gphmmdim(){

   return phmm2ph.size();

}
