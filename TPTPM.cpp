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

int **TPTPM::t2tpmm;
vector< vector<int> > TPTPM::tpmm2t;

/**
 * initialize the static lists
 */
void TPTPM::init(){

   t2tpmm = new int * [TPM::gn()];

   for(int i = 0;i < TPM::gn();++i)
      t2tpmm[i] = new int [TPM::gn()];

   vector<int> v(2);

   int tpmm = 0;

   for(int i = 0;i < TPM::gn();++i)
      for(int j = i;j < TPM::gn();++j){

         v[0] = i;
         v[1] = j;

         tpmm2t.push_back(v);

         t2tpmm[i][j] = tpmm;
         t2tpmm[j][i] = tpmm;

         ++tpmm;

      }

}

/**
 * deallocate the static lists
 */
void TPTPM::clear(){

   for(int i = 0;i < TPM::gn();++i)
      delete [] t2tpmm[i];

   delete [] t2tpmm;

}

/**
 * standard constructor:
 */
TPTPM::TPTPM() : Matrix(tpmm2t.size()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix tpmm_c
 * @param tpmm_c object that will be copied into this.
 */
TPTPM::TPTPM(const TPTPM &tpmm_c) : Matrix(tpmm_c){ }

/**
 * destructor
 */
TPTPM::~TPTPM(){ }

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
 * @return the number on place TPTPM(i,j) with the right phase.
 */
double TPTPM::operator()(int a,int b,int c,int d,int e,int z,int t,int h) const{

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
 * access the elements of the matrix in tp mode
 * @param I first tp index that forms the tpmm row index i together with J
 * @param J second tp index that forms the tpmm row index i together with I
 * @param K first tp index that forms the tpmm column index j together with L
 * @param L second tp index that forms the tpmm column index j together with K
 * @return the number on place TPTPM(i,j)
 */
double TPTPM::operator()(int I,int J,int K,int L) const{

   int i = t2tpmm[I][J];
   int j = t2tpmm[K][L];

   return (*this)(i,j);

}

ostream &operator<<(ostream &output,const TPTPM &tpmm_p){

   int I,J,K,L;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      I = tpmm_p.tpmm2t[i][0];
      J = tpmm_p.tpmm2t[i][1];

      a = TPM::gt2s(I,0);
      b = TPM::gt2s(I,1);

      c = TPM::gt2s(J,0);
      d = TPM::gt2s(J,1);

      for(int j = i;j < TPTPM::gn();++j){

         K = tpmm_p.tpmm2t[j][0];
         L = tpmm_p.tpmm2t[j][1];

         e = TPM::gt2s(K,0);
         z = TPM::gt2s(K,1);

         t = TPM::gt2s(L,0);
         h = TPM::gt2s(L,1);

         output << i << "\t" << j << "\t|\t" << I << "\t" << J << "\t" << K << "\t" << L << "\t|\t" << 
         
            "(" << a << "," << b << "," << c << "," << d << ")\t(" << e << "," << z << "," << t << "," << h << ")\t|\t" << tpmm_p(i,j) << endl;

      }

   }

   return output;

}

/**
 * @return the dimension of a TPTPM matrix
 */
int TPTPM::gn(){

   return tpmm2t.size();

}

/**
 * access to the lists from outside the class
 */
int TPTPM::gt2tpmm(int i,int j){

   return t2tpmm[i][j];

}

/**
 * access to the lists from outside the class
 * @param option == 0 return a, == 1 return b
 */
int TPTPM::gtpmm2t(int i,int option){

   return tpmm2t[i][option];

}

 /**
 * construct the antisymmetrized "symmetric direct product" of two PHM matrix
 */
void TPTPM::dp(const PHM &phm){ 

   int I,J,K,L;

   int a,b,c,d,e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      I = TPTPM::gtpmm2t(i,0);
      J = TPTPM::gtpmm2t(i,1);

      a = TPM::gt2s(I,0);
      b = TPM::gt2s(I,1);

      c = TPM::gt2s(J,0);
      d = TPM::gt2s(J,1);

      for(int j = i;j < TPTPM::gn();++j){

         K = TPTPM::gtpmm2t(j,0);
         L = TPTPM::gtpmm2t(j,1);

         e = TPM::gt2s(K,0);
         z = TPM::gt2s(K,1);

         t = TPM::gt2s(L,0);
         h = TPM::gt2s(L,1);

         //first 16 direct product terms of the PHM object
         (*this)(i,j) += phm(a,d,e,h) * phm(c,b,t,z) + phm(a,d,t,z) * phm(c,b,e,h) - phm(a,d,z,h) * phm(c,b,t,e) - phm(a,d,t,e) * phm(c,b,z,h)

               - phm(a,d,e,t) * phm(c,b,h,z) - phm(a,d,h,z) * phm(c,b,e,t) + phm(a,d,z,t) * phm(c,b,h,e) + phm(a,d,h,e) * phm(c,b,z,t)

               - phm(b,d,e,h) * phm(c,a,t,z) - phm(b,d,t,z) * phm(c,a,e,h) + phm(b,d,z,h) * phm(c,a,t,e) + phm(b,d,t,e) * phm(c,a,z,h)

               + phm(b,d,e,t) * phm(c,a,h,z) + phm(b,d,h,z) * phm(c,a,e,t) - phm(b,d,z,t) * phm(c,a,h,e) - phm(b,d,h,e) * phm(c,a,z,t)

               - phm(a,c,e,h) * phm(d,b,t,z) - phm(a,c,t,z) * phm(d,b,e,h) + phm(a,c,z,h) * phm(d,b,t,e) + phm(a,c,t,e) * phm(d,b,z,h)

               + phm(a,c,e,t) * phm(d,b,h,z) + phm(a,c,h,z) * phm(d,b,e,t) - phm(a,c,z,t) * phm(d,b,h,e) - phm(a,c,h,e) * phm(d,b,z,t)

               + phm(b,c,e,h) * phm(d,a,t,z) + phm(b,c,t,z) * phm(d,a,e,h) - phm(b,c,z,h) * phm(d,a,t,e) - phm(b,c,t,e) * phm(d,a,z,h)

               - phm(b,c,e,t) * phm(d,a,h,z) - phm(b,c,h,z) * phm(d,a,e,t) + phm(b,c,z,t) * phm(d,a,h,e) + phm(b,c,h,e) * phm(d,a,z,t) ;

      }

   }

   this->symmetrize();

}

/**
 * construct the antisymmetrized, double 'tilde' of a symmetric direct product of two PPHM matrices
 * @param pphm input PPHM
 */
void TPTPM::dpw2_slow(const PPHM &pphm){ 

   int M = Tools::gM();

   int I,J,K,L;

   int a,b,c,d,e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      I = TPTPM::gtpmm2t(i,0);
      J = TPTPM::gtpmm2t(i,1);

      a = TPM::gt2s(I,0);
      b = TPM::gt2s(I,1);

      c = TPM::gt2s(J,0);
      d = TPM::gt2s(J,1);

      for(int j = i;j < TPTPM::gn();++j){

         K = TPTPM::gtpmm2t(j,0);
         L = TPTPM::gtpmm2t(j,1);

         e = TPM::gt2s(K,0);
         z = TPM::gt2s(K,1);

         t = TPM::gt2s(L,0);
         h = TPM::gt2s(L,1);

         (*this)(i,j) = 0.0;

         for(int k = 0;k < M;++k)
            for(int l = 0;l < M;++l){

               (*this)(i,j) += pphm(k,d,a,l,h,e) * pphm(k,b,c,l,z,t) + pphm(k,d,a,l,z,t) * pphm(k,b,c,l,h,e)

                  - pphm(k,d,b,l,h,e) * pphm(k,a,c,l,z,t) - pphm(k,d,b,l,z,t) * pphm(k,a,c,l,h,e)

                  - pphm(k,c,a,l,h,e) * pphm(k,b,d,l,z,t) - pphm(k,c,a,l,z,t) * pphm(k,b,d,l,h,e)

                  + pphm(k,c,b,l,h,e) * pphm(k,a,d,l,z,t) + pphm(k,c,b,l,z,t) * pphm(k,a,d,l,h,e)
 
                  - pphm(k,d,a,l,h,z) * pphm(k,b,c,l,e,t) - pphm(k,d,a,l,e,t) * pphm(k,b,c,l,h,z)

                  + pphm(k,d,b,l,h,z) * pphm(k,a,c,l,e,t) + pphm(k,d,b,l,e,t) * pphm(k,a,c,l,h,z)

                  + pphm(k,c,a,l,h,z) * pphm(k,b,d,l,e,t) + pphm(k,c,a,l,e,t) * pphm(k,b,d,l,h,z)

                  - pphm(k,c,b,l,h,z) * pphm(k,a,d,l,e,t) - pphm(k,c,b,l,e,t) * pphm(k,a,d,l,h,z)

                  - pphm(k,d,a,l,t,e) * pphm(k,b,c,l,z,h) - pphm(k,d,a,l,z,h) * pphm(k,b,c,l,t,e)

                  + pphm(k,d,b,l,t,e) * pphm(k,a,c,l,z,h) + pphm(k,d,b,l,z,h) * pphm(k,a,c,l,t,e)

                  + pphm(k,c,a,l,t,e) * pphm(k,b,d,l,z,h) + pphm(k,c,a,l,z,h) * pphm(k,b,d,l,t,e)

                  - pphm(k,c,b,l,t,e) * pphm(k,a,d,l,z,h) - pphm(k,c,b,l,z,h) * pphm(k,a,d,l,t,e)

                  + pphm(k,d,a,l,t,z) * pphm(k,b,c,l,e,h) + pphm(k,d,a,l,e,h) * pphm(k,b,c,l,t,z)

                  - pphm(k,d,b,l,t,z) * pphm(k,a,c,l,e,h) - pphm(k,d,b,l,e,h) * pphm(k,a,c,l,t,z)

                  - pphm(k,c,a,l,t,z) * pphm(k,b,d,l,e,h) - pphm(k,c,a,l,e,h) * pphm(k,b,d,l,t,z)

                  + pphm(k,c,b,l,t,z) * pphm(k,a,d,l,e,h) + pphm(k,c,b,l,e,h) * pphm(k,a,d,l,t,z);
 
            }

      }
   }

   this->symmetrize();

}

/**
 * yet another double trace of a direct product of DPM's
 */
void TPTPM::dpt2(const DPM &dpm){

   int M = Tools::gM();
   int M2 = M*M;
   int M3 = M2*M;
   int M4 = M3*M;
   int M5 = M4*M;

   double *dparray = new double [M5*M];

   dpm.convert(dparray);

   int I,J,K,L;

   int a,b,c,d,e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      I = TPTPM::gtpmm2t(i,0);
      J = TPTPM::gtpmm2t(i,1);

      a = TPM::gt2s(I,0);
      b = TPM::gt2s(I,1);

      c = TPM::gt2s(J,0);
      d = TPM::gt2s(J,1);

      for(int j = i;j < TPTPM::gn();++j){

         K = TPTPM::gtpmm2t(j,0);
         L = TPTPM::gtpmm2t(j,1);

         e = TPM::gt2s(K,0);
         z = TPM::gt2s(K,1);

         t = TPM::gt2s(L,0);
         h = TPM::gt2s(L,1);

         (*this)(i,j) = 0.0;

         for(int k = 0;k < M;++k)
            for(int l = 0;l < M;++l){

               (*this)(i,j) += dparray[a*M5 + b*M4 + k*M3 + e*M2 + z*M + l] * dparray[c*M5 + d*M4 + k*M3 + t*M2 + h*M + l]

                  + dparray[a*M5 + b*M4 + k*M3 + t*M2 + h*M + l] * dparray[c*M5 + d*M4 + k*M3 + e*M2 + z*M + l]; 

            }

      }
   }

   this->symmetrize();

   delete [] dparray;

}
