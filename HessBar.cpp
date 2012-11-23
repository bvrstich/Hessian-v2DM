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

int ***HessBar::s2hb;
vector< vector<int> > HessBar::hb2s;

/**
 * initialize the static lists
 */
void HessBar::init(){

   s2hb = new int ** [Tools::gM()];

   for(int a = 0;a < Tools::gM();++a){

      s2hb[a] = new int * [Tools::gM()];

      for(int b = 0;b < Tools::gM();++b)
         s2hb[a][b] = new int [Tools::gM()];

   }

   vector<int> v(3);

   int hb = 0;

   for(int a = 0;a < Tools::gM();++a)
      for(int b = a + 1;b < Tools::gM();++b)
         for(int c = 0;c < Tools::gM();++c){

            v[0] = a;
            v[1] = b;
            v[2] = c;

            hb2s.push_back(v);

            s2hb[a][b][c] = hb;
            s2hb[b][a][c] = hb;

            ++hb;

         }

}

/**
 * deallocate the static lists
 */
void HessBar::clear(){

   for(int a = 0;a < Tools::gM();++a){

      for(int b = 0;b < Tools::gM();++b)
         delete [] s2hb[a][b];

      delete [] s2hb[a];

   }

   delete [] s2hb;

}

/**
 * standard constructor:
 */
HessBar::HessBar() : Matrix(hb2s.size()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix hb_c
 * @param hb_c object that will be copied into this.
 */
HessBar::HessBar(const HessBar &hb_c) : Matrix(hb_c){ }

/**
 * destructor
 */
HessBar::~HessBar(){ }

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:\n\n
 * @param a first sp index that forms the tp row index i together with b
 * @param b second sp index that forms the tp row index i together with a
 * @param c first sp index that forms the tp column index j together with d
 * @param d second sp index that forms the tp column index j together with c
 * @param e first sp index that forms the tp row index i together with b
 * @param z second sp index that forms the tp row index i together with a
 * @return the number on place HessBar(i,j) with the right phase.
 */
double HessBar::operator()(int a,int b,int c,int d,int e,int z) const{

   if( (a == b) || (d == e) )
      return 0;
   else{

      int phase = 1;

      if(a > b)
         phase *= -1;

      if(d > e)
         phase *= -1;

      int i = s2hb[a][b][c];
      int j = s2hb[d][e][z];

      return phase * (*this)(i,j);

   }

}

ostream &operator<<(ostream &output,const HessBar &hb_p){

   for(int i = 0;i < hb_p.gn();++i)
      for(int j = 0;j < hb_p.gn();++j){

         output << i << "\t" << j << "\t|\t" << hb_p.hb2s[i][0] << "\t" << hb_p.hb2s[i][1] << "\t" << hb_p.hb2s[i][2]

            << "\t" << hb_p.hb2s[j][0] << "\t" << hb_p.hb2s[j][1] << "\t" << hb_p.hb2s[j][2] << "\t" << hb_p(i,j) << endl;

      }

   return output;

}

/**
 * access to the lists from outside the class
 */
int HessBar::gs2hb(int a,int b,int c){

   return s2hb[a][b][c];

}

/**
 * access to the lists from outside the class
 * @param option == 0 return a, == 1 return b, == 2 return c
 */
int HessBar::ghb2s(int i,int option){

   return hb2s[i][option];

}

/**
 * fill the HessBar object by tracing out one pair of indices of the direct product of a TPM object with itsself
 * @param scale factor to scale the trace with
 */
void HessBar::dirprodtrace(double scale,const TPM &Q){

   int a,b,c,d,e,z;

   for(int i = 0;i < gn();++i){

      a = hb2s[i][0];
      b = hb2s[i][1];
      c = hb2s[i][2];

      for(int j = i;j < gn();++j){

         d = hb2s[j][0];
         e = hb2s[j][1];
         z = hb2s[j][2];

         (*this)(i,j) = 0.0;

         for(int l = 0;l < Tools::gM();++l)
            (*this)(i,j) += Q(a,b,c,l) * Q(d,e,z,l) + Q(d,e,c,l) * Q(a,b,z,l);

         (*this)(i,j) *= scale;

      }
   }

   this->symmetrize();

}
