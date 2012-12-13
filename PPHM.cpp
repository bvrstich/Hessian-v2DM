#include <iostream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

vector< vector<int> > PPHM::pph2s;
int ***PPHM::s2pph;

/**
 * initialize the static variables and allocate the static lists
 */
void PPHM::init(){

   int M = Tools::gM();

   s2pph = new int ** [M];

   for(int a = 0;a < M;++a){

      s2pph[a] = new int * [M];

      for(int b = 0;b < M;++b)
         s2pph[a][b] = new int [M];

   }

   int pph = 0;

   vector<int> v(3);

   for(int a = 0;a < M;++a)
      for(int b = a + 1;b < M;++b)
         for(int c = 0;c < M;++c){

            s2pph[a][b][c] = pph;

            v[0] = a;
            v[1] = b;
            v[2] = c;

            pph2s.push_back(v);

            ++pph;

         }

}

/**
 * deallocate the static lists
 */
void PPHM::clear(){

   for(int a = 0;a < Tools::gM();++a){

      for(int b = 0;b < Tools::gM();++b)
         delete [] s2pph[a][b];

      delete [] s2pph[a];

   }

   delete [] s2pph;

}

/**
 * standard constructor: constructs Matrix object of dimension M*M*(M - 1)/2 and
 */
PPHM::PPHM() : Matrix(pph2s.size()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*M*(M - 1)/2 and copies the content of pphm_c into it,
 * @param pphm_c input PPHM to be copied
 */
PPHM::PPHM(const PPHM &pphm_c) : Matrix(pphm_c){ }

/**
 * destructor
 */
PPHM::~PPHM(){ }

/**
 * access the elements of the matrix in sp mode, antisymmetry in the first two indices is automatically accounted for:\n\n
 * PPHM(a,b,c,d,e,f) = -PPHM(b,a,c,d,e,f)  but PPHM(a,b,c,d,e,f) != - PPHM(a,c,b,d,e,f) \n\n
 * PPHM(a,a,c,d,e,f) = 0 but PPHM(a,b,b,d,e,f) != 0 \n\n
 * @param a first sp index that forms the pph row index i together with b and c
 * @param b second sp index that forms the pph row index i together with a and c
 * @param c third sp index that forms the pph row index i together with a and b
 * @param d first sp index that forms the pph column index j together with e and z
 * @param e second sp index that forms the pph column index j together with d and z
 * @param z third sp index that forms the pph column index j together with d and e
 * @return the number on place PPHM(i,j) with the right phase.
 */
double PPHM::operator()(int a,int b,int c,int d,int e,int z) const{

   //eerst kijken of de eerste twee indices gelijk zijn:
   if(a == b || d == e)
      return 0;

   //dan kijken wel pph index met welke fase moet genomen worden:
   //eerst voor de i
   int i;

   int phase = 1;

   if(a > b){

      i = s2pph[b][a][c];
      phase *= -1;

   }
   else
      i = s2pph[a][b][c];

   int j;

   if(d > e){

      j = s2pph[e][d][z];
      phase *= -1;

   }
   else
      j = s2pph[d][e][z];

   return phase*(*this)(i,j);

}

ostream &operator<<(ostream &output,PPHM &pphm_p){

   for(int i = 0;i < pphm_p.gn();++i)
      for(int j = 0;j < pphm_p.gn();++j){

         output << i << "\t" << j << "\t|\t" << pphm_p.pph2s[i][0] << "\t" << pphm_p.pph2s[i][1] << "\t" << pphm_p.pph2s[i][2]

            << "\t" << pphm_p.pph2s[j][0] << "\t" << pphm_p.pph2s[j][1] << "\t" << pphm_p.pph2s[j][2] << "\t" << pphm_p(i,j) << endl;

      }

   return output;

}

/**
 * The T2-map: maps a TPM object (tpm) on a PPHM object (*this). see primal_dual.pdf for more information
 * @param option == 0, regular T2, == 1, special (incorrect) T2, keep for test in program with regular T2
 * @param tpm input TPM
 */
void PPHM::T(int option,const TPM &tpm){

   if(option == 0){

      //construct the spm
      SPM spm;
      spm.bar(1.0/(Tools::gN() - 1.0),tpm);

      int a,b,c,d,e,z;

      for(int i = 0;i < gn();++i){

         a = pph2s[i][0];
         b = pph2s[i][1];
         c = pph2s[i][2];

         for(int j = i;j < gn();++j){

            d = pph2s[j][0];
            e = pph2s[j][1];
            z = pph2s[j][2];

            //initialize
            (*this)(i,j) = 0.0;

            if(a == d){

               //sp part
               if(b == e)
                  (*this)(i,j) += spm(c,z);

               //tp part
               (*this)(i,j) -= tpm(c,e,z,b);

            }

            //now only tp parts left:
            if(c == z)
               (*this)(i,j) += tpm(a,b,d,e);

            if(b == d)
               (*this)(i,j) += tpm(c,e,z,a);

            if(b == e)
               (*this)(i,j) -= tpm(c,d,z,a);

         }
      }
   }
   else{

      TPM Q;
      Q.Q(1,tpm);

      DPM T;
      T.T(1,tpm);

      int a,b,c,d,e,z;

      for(int i = 0;i < gn();++i){

         a = pph2s[i][0];
         b = pph2s[i][1];
         c = pph2s[i][2];

         for(int j = i;j < gn();++j){

            d = pph2s[j][0];
            e = pph2s[j][1];
            z = pph2s[j][2];

            //initialize
            (*this)(i,j) = 0.0;

            if(c == z)
               (*this)(i,j) += tpm(a,b,d,e) + Q(a,b,d,e);

            (*this)(i,j) -= T(a,b,z,d,e,c);

         }
      }

   }

   //and symmetrize
   this->symmetrize();

}

/**
 * access the numbers of the PPHM matrix in sp mode, but where the input sp indices
 * are ordered (a < b < c, d < e < z) so that no if statements and other bullcrap is needed.
 */
const double &PPHM::ordacc(int a,int b,int c,int d,int e,int z) const {
   
   return  (*this)(s2pph[a][b][c],s2pph[d][e][z]);

}

/** 
 * yet another way to avoid the very slow access operator with the ordering
 */
void PPHM::convert(double *ppharray) const{

   int M = Tools::gM();
   int M2 = M*M;
   int M3 = M2*M;
   int M4 = M3*M;
   int M5 = M4*M;

   for(int a = 0;a < M;++a)
      for(int b = 0;b < M;++b)
         for(int c = 0;c < M;++c)
            for(int d = 0;d < M;++d)
               for(int e = 0;e < M;++e)
                  for(int z = 0;z < M;++z)
                     ppharray[a*M5 + b*M4 + c*M3 + d*M2 + e*M + z] = (*this)(a,b,c,d,e,z);


}

/* vim: set ts=3 sw=3 expandtab :*/
