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

int **TPTPnsM::t2tpmm;
vector< vector<int> > TPTPnsM::tpmm2t;

/**
 * initialize the static lists
 */
void TPTPnsM::init(){

   t2tpmm = new int * [TPM::gn()];

   for(int i = 0;i < TPM::gn();++i)
      t2tpmm[i] = new int [TPM::gn()];

   vector<int> v(2);

   int tpmm = 0;

   //no symmetry, difference with t2hess
   for(int i = 0;i < TPM::gn();++i)
      for(int j = 0;j < TPM::gn();++j){

         v[0] = i;
         v[1] = j;

         tpmm2t.push_back(v);

         t2tpmm[i][j] = tpmm;

         ++tpmm;

      }

}

/**
 * deallocate the static lists
 */
void TPTPnsM::clear(){

   for(int i = 0;i < TPM::gn();++i)
      delete [] t2tpmm[i];

   delete [] t2tpmm;

}

/**
 * standard constructor:
 */
TPTPnsM::TPTPnsM() : Matrix(tpmm2t.size()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix tpmm_c
 * @param tpmm_c object that will be copied into this.
 */
TPTPnsM::TPTPnsM(const TPTPnsM &tpmm_c) : Matrix(tpmm_c){ }

/**
 * destructor
 */
TPTPnsM::~TPTPnsM(){ }

/**
 * access the elements of the matrix in tp mode
 * @param I first tp index that forms the tpmm row index i together with J
 * @param J second tp index that forms the tpmm row index i together with I
 * @param K first tp index that forms the tpmm column index j together with L
 * @param L second tp index that forms the tpmm column index j together with K
 * @return the number on place TPTPnsM(i,j) with the right phase.
 */
double TPTPnsM::operator()(int I,int J,int K,int L) const {

   int i = t2tpmm[I][J];
   int j = t2tpmm[K][L];

   return (*this)(i,j);

}

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:\n\n
 * @param a first sp index that forms the tp row index i together with b
 * @param b second sp index that forms the tp row index i together with a
 * @param c first sp index that forms the tp column index j together with d
 * @param d second sp index that forms the tp column index j together with c
 * @param e first sp index that forms the tp row index i together with b
 * @param z second sp index that forms the tp row index i together with a
 * @param t first sp index that forms the tp column index j together with d
 * @param h second sp index that forms the tp column index j together with c
 * @return the number on place Hessian(i,j) with the right phase.
 */
double TPTPnsM::operator()(int a,int b,int c,int d,int e,int z,int t,int h) const {

   if( (a == b) || (c == d) || (e == z) || (t == h))
      return 0;
   else{

      int I = TPM::gs2t(a,b);
      int J = TPM::gs2t(c,d);
      int K = TPM::gs2t(e,z);
      int L = TPM::gs2t(t,h);

      int phase = 1;

      if(a > b)
         phase *= -1;

      if(c > d)
         phase *= -1;

      if(e > z)
         phase *= -1;

      if(t > h)
         phase *= -1;

      int i = t2tpmm[I][J];
      int j = t2tpmm[K][L];

      return phase * (*this)(i,j);

   }

}

/**
 * square a TPSPM object to obtain a TPTPnsM
 * @param tpspm input TPSPM
 */
void TPTPnsM::square(const TPSPnsM &tpspnsm){

   char transA = 'N';
   char transB = 'T';

   double alpha = 1.0;
   double beta = 0.0;

   int m = tpspnsm.gm();
   int n = tpspnsm.gn();

   dgemm_(&transA,&transB,&m,&m,&n,&alpha,tpspnsm.gRecMat()[0],&m,tpspnsm.gRecMat()[0],&m,&beta,this->gMatrix()[0],&m);

}

/**
 * @return the dimension of the matrix
 */
int TPTPnsM::gn(){

   return tpmm2t.size();

}

/**
 * access to the lists from outside of the class
 */
int TPTPnsM::gtpmm2t(int i,int option){

   return tpmm2t[i][option];

}

/**
 * access to the lists from outside the class
 */
int TPTPnsM::gt2tpmm(int i,int j){

   return t2tpmm[i][j];

}
