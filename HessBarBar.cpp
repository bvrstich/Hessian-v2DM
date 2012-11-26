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
HessBarBar::HessBarBar() : Matrix(HessBar::gncol()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix hbb_c
 * @param hbb_c object that will be copied into this.
 */
HessBarBar::HessBarBar(const HessBarBar &hbb_c) : Matrix(hbb_c){ }

/**
 * destructor
 */
HessBarBar::~HessBarBar(){ }

/**
 * access the elements of the matrix in sp mode
 * @param a first sp index that forms the hbb row index i together with b
 * @param b second sp index that forms the hbb row index i together with a
 * @param c first sp index that forms the hbb column index j together with d
 * @param d second sp index that forms the hbb column index j together with c
 * @return the number on place HessBarBar(i,j) with the right phase.
 */
double HessBarBar::operator()(int a,int b,int c,int d) const {

   int i = HessBar::gs2hb(a,b);
   int j = HessBar::gs2hb(c,d);

   return (*this)(i,j);

}

/**
 * construct a HessBarBar object by tracing one pair of indices of a HessBar object
 * @param scale scalefactor for the barred object
 * @param hb input HessBar object
 */
void HessBarBar::bar(double scale,const HessBar &hb){

   int a,b,c,d;

   for(int i = 0;i < gn();++i){

      a = HessBar::ghb2s(i,0);
      b = HessBar::ghb2s(i,1);

      for(int j = i;j < gn();++j){

         c = HessBar::ghb2s(j,0);
         d = HessBar::ghb2s(j,1);

         (*this)(i,j) = 0.0;

         for(int l = 0;l < Tools::gM();++l)
            (*this)(i,j) += hb(a,l,b,l,c,d);

         (*this)(i,j) *= scale;

      }
   }

   this->symmetrize();

}
