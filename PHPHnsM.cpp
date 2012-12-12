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
PHPHnsM::PHPHnsM() : Matrix(PHSPnsM::gphmmdim()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix phmm_c
 * @param phmm_c object that will be copied into this.
 */
PHPHnsM::PHPHnsM(const PHPHnsM &phmm_c) : Matrix(phmm_c){ }

/**
 * destructor
 */
PHPHnsM::~PHPHnsM(){ }

/**
 * access the elements of the matrix in ph mode
 * @param I first ph index that forms the phmm row index i together with J
 * @param J second ph index that forms the phmm row index i together with I
 * @param K first ph index that forms the phmm column index j together with L
 * @param L second ph index that forms the phmm column index j together with K
 * @return the number on place PHPHnsM(i,j) with the right phase.
 */
double PHPHnsM::operator()(int I,int J,int K,int L) const {

   int i = PHSPnsM::gph2phmm(I,J);
   int j = PHSPnsM::gph2phmm(K,L);

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
 * @return the number on place PHPHnsM(i,j)
 */
double PHPHnsM::operator()(int a,int b,int c,int d,int e,int z,int t,int h) const {

   int I = PHM::gs2ph(a,b);
   int J = PHM::gs2ph(c,d);
   int K = PHM::gs2ph(e,z);
   int L = PHM::gs2ph(t,h);

   int i = PHSPnsM::gph2phmm(I,J);
   int j = PHSPnsM::gph2phmm(K,L);

   return (*this)(i,j);

}

/**
 * square a PHSPnsM object to obtain a PHPHnsM
 * @param phspnsm input PHSPnsM
 */
void PHPHnsM::square(const PHSPnsM &phspnsm){

   char transA = 'N';
   char transB = 'T';

   double alpha = 1.0;
   double beta = 0.0;

   int m = phspnsm.gm();
   int n = phspnsm.gn();

   cout << m << "\t" << n << endl;

   dgemm_(&transA,&transB,&m,&m,&n,&alpha,phspnsm.gRecMat()[0],&m,phspnsm.gRecMat()[0],&m,&beta,this->gMatrix()[0],&m);

}
