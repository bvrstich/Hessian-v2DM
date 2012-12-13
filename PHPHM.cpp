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

vector< vector<int> > PHPHM::phmm2ph;
int **PHPHM::ph2phmm;

/**
 * initialize the static lists
 */
void PHPHM::init(){

  ph2phmm = new int * [PHM::gn()];

   for(int i = 0;i < PHM::gn();++i)
      ph2phmm[i] = new int [PHM::gn()];

   vector<int> v(2);

   int phmm = 0;

   //symmetry between i and j!
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
void PHPHM::clear(){

   for(int i = 0;i < PHM::gn();++i)
      delete [] ph2phmm[i];

   delete [] ph2phmm;

}

/**
 * standard constructor:
 */
PHPHM::PHPHM() : Matrix(phmm2ph.size()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix phmm_c
 * @param phmm_c object that will be copied into this.
 */
PHPHM::PHPHM(const PHPHM &phmm_c) : Matrix(phmm_c){ }

/**
 * destructor
 */
PHPHM::~PHPHM(){ }

/**
 * access the elements of the matrix in ph mode
 * @param I first ph index that forms the phmm row index i together with J
 * @param J second ph index that forms the phmm row index i together with I
 * @param K first ph index that forms the phmm column index j together with L
 * @param L second ph index that forms the phmm column index j together with K
 * @return the number on place PHPHM(i,j) with the right phase.
 */
double PHPHM::operator()(int I,int J,int K,int L) const {

   int i = ph2phmm[I][J];
   int j = ph2phmm[K][L];

   return (*this)(i,j);

}

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:\n\n
 * @param a first sp index that forms the ph row index i together with b
 * @param b second sp index that forms the ph row index i together with a
 * @param c first sp index that forms the ph column index j together with d
 * @param d second sp index that forms the ph column index j together with c
 * @param e first sp index that forms the ph row index i together with b
 * @param z second sp index that forms the ph row index i together with a
 * @param t first sp index that forms the ph column index j together with d
 * @param h second sp index that forms the ph column index j together with c
 * @return the number on place PHPHM(i,j)
 */
double PHPHM::operator()(int a,int b,int c,int d,int e,int z,int t,int h) const {

   int I = PHM::gs2ph(a,b);
   int J = PHM::gs2ph(c,d);
   int K = PHM::gs2ph(e,z);
   int L = PHM::gs2ph(t,h);

   int i = ph2phmm[I][J];
   int j = ph2phmm[K][L];

   return (*this)(i,j);

}
