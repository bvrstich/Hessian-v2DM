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
SPSPM::SPSPM() : Matrix(TPSPM::gncol()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix spspm_c
 * @param spspm_c object that will be copied into this.
 */
SPSPM::SPSPM(const SPSPM &spspm_c) : Matrix(spspm_c){ }

/**
 * destructor
 */
SPSPM::~SPSPM(){ }

/**
 * access the elements of the matrix in sp mode
 * @param a first sp index that forms the spspm row index i together with b
 * @param b second sp index that forms the spspm row index i together with a
 * @param c first sp index that forms the spspm column index j together with d
 * @param d second sp index that forms the spspm column index j together with c
 * @return the number on place SPSPM(i,j) with the right phase.
 */
double SPSPM::operator()(int a,int b,int c,int d) const {

   int i = TPSPM::gs2spmm(a,b);
   int j = TPSPM::gs2spmm(c,d);

   return (*this)(i,j);

}

/**
 * construct a SPSPM object by tracing one pair of indices of a TPSPM object
 * @param scale scalefactor for the barred object
 * @param hb input TPSPM object
 */
void SPSPM::bar(double scale,const TPSPM &hb){

   int a,b,c,d;

   for(int i = 0;i < gn();++i){

      a = TPSPM::gspmm2s(i,0);
      b = TPSPM::gspmm2s(i,1);

      for(int j = i;j < gn();++j){

         c = TPSPM::gspmm2s(j,0);
         d = TPSPM::gspmm2s(j,1);

         (*this)(i,j) = 0.0;

         for(int l = 0;l < Tools::gM();++l)
            (*this)(i,j) += hb(a,l,b,l,c,d);

         (*this)(i,j) *= scale;

      }
   }

   this->symmetrize();

}
