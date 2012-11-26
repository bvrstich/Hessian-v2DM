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

int **HessBar::s2hb;
vector< vector<int> > HessBar::hb2s;

/**
 * initialize the static lists
 */
void HessBar::init(){

   s2hb = new int * [Tools::gM()];

   for(int a = 0;a < Tools::gM();++a)
      s2hb[a] = new int [Tools::gM()];

   vector<int> v(2);

   int hb = 0;

   for(int a = 0;a < Tools::gM();++a)
      for(int b = a;b < Tools::gM();++b){

         v[0] = a;
         v[1] = b;

         hb2s.push_back(v);

         s2hb[a][b] = hb;
         s2hb[b][a] = hb;

         ++hb;

      }

}

/**
 * deallocate the static lists
 */
void HessBar::clear(){

   for(int a = 0;a < Tools::gM();++a)
      delete [] s2hb[a];

   delete [] s2hb;

}

/**
 * standard constructor:
 */
HessBar::HessBar() : RecMat(Hessian::gn(),hb2s.size()) { }

/**
 * copy constructor: constructs RecMat object of dimension M*(M - 1)/2 and fills it with the content of matrix hb_c
 * @param hb_c object that will be copied into this.
 */
HessBar::HessBar(const HessBar &hb_c) : RecMat(hb_c){ }

/**
 * destructor
 */
HessBar::~HessBar(){ }

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:\n\n
 * @param a first sp index that forms the tp row index i together with b, part of the HessBar row index
 * @param b second sp index that forms the tp row index i together with a, part of the HessBar row index
 * @param c first sp index that forms the tp column index j together with d, part of the HessBar row index
 * @param d second sp index that forms the tp column index j together with c, part of the HessBar row index
 * @param e first sp index that forms part of the HessBar colum index together with b
 * @param z second sp index that forms part of the HessBar colum index together with a
 * @return the number on place HessBar(i,j) with the right phase.
 */
double HessBar::operator()(int a,int b,int c,int d,int e,int z) const{

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

      int i = Hessian::gt2hess(I,J);

      int j = s2hb[e][z];

      return phase * (*this)(i,j);

   }

}

/**
 * access to the column index list from outside the class
 */
int HessBar::gs2hb(int a,int b){

   return s2hb[a][b];

}

/**
 * access to the column index list from outside the class
 * @param option == 0 return a, == 1 return b
 */
int HessBar::ghb2s(int i,int option){

   return hb2s[i][option];

}

/**
 * fill the HessBar object by tracing out one pair of indices of the direct product of a TPM object with itsself
 * @param scale factor to scale the trace with
 */
void HessBar::dirprodtrace(double scale,const TPM &Q){

   int I,J;

   int a,b,c,d,e,z;

   for(int i = 0;i < gm();++i){

      I = Hessian::ghess2t(i,0);
      J = Hessian::ghess2t(i,1);

      a = TPM::gt2s(I,0);
      b = TPM::gt2s(I,1);

      c = TPM::gt2s(J,0);
      d = TPM::gt2s(J,1);

      for(int j = 0;j < gn();++j){

         e = hb2s[j][0];
         z = hb2s[j][1];

         (*this)(i,j) = 0.0;

         for(int l = 0;l < Tools::gM();++l)
            (*this)(i,j) += Q(a,b,e,l) * Q(c,d,z,l) + Q(c,d,e,l) * Q(a,b,z,l);

         (*this)(i,j) *= scale;

      }
   }

}

/**
 * @return the number of columns
 */
int HessBar::gncol(){

   return hb2s.size();

}
