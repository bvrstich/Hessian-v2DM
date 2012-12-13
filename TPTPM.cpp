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
 * Fills a TPTPM object with the two-times traced symmetric direct product of a DPM object.
 * @param dpm input DPM object
 */
void TPTPM::dirprodtrace(const DPM &dpm){

   int M = Tools::gM();

   int I,J,K,L;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < gn();++i){

      I = tpmm2t[i][0];
      J = tpmm2t[i][1];

      a = TPM::gt2s(I,0);
      b = TPM::gt2s(I,1);

      c = TPM::gt2s(J,0);
      d = TPM::gt2s(J,1);

      if(b <= c){

         for(int j = i;j < gn();++j){

            K = tpmm2t[j][0];
            L = tpmm2t[j][1];

            e = TPM::gt2s(K,0);
            z = TPM::gt2s(K,1);

            t = TPM::gt2s(L,0);
            h = TPM::gt2s(L,1);

            (*this)(i,j) = 0.0;

            if(z <= t){

               for(int k = 0;k < a;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(k,a,b,t,l,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(k,a,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = a + 1;k < b;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = b + 1;k < c;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = c + 1;k < d;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,k,d,l,e,z);

                  for(int l = e + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = z + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,k,d,e,z,l);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,k,d,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,k,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,k,d,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,k,d,e,z,l);

               }

               for(int k = d + 1;k < M;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,l,e,z);

                  for(int l = e + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = z + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,e,z,l);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,d,k,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,d,k,e,z,l);

               }

            }
            else if(z <= h){

               for(int k = 0;k < a;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = t + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(k,a,b,t,l,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(k,a,b,t,l,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(k,a,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = a + 1;k < b;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = t + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = b + 1;k < c;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = t + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = c + 1;k < d;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,k,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = t + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,k,d,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = z + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,k,d,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,k,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,k,d,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,k,d,e,z,l);

               }

               for(int k = d + 1;k < M;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = t + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = z + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,d,k,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,d,k,e,z,l);

               }

            }
            else{

               for(int k = 0;k < a;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(k,a,b,t,l,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = h + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(k,a,b,t,h,l) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(k,a,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = a + 1;k < b;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = h + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = b + 1;k < c;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = h + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = c + 1;k < d;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,k,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,k,d,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = h + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,k,d,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = z + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,k,d,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,k,d,e,z,l);

               }

               for(int k = d + 1;k < M;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = h + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = z + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,d,k,e,z,l);

               }

            }

         }

      }
      else if(b <= d){

         for(int j = i;j < gn();++j){

            K = tpmm2t[j][0];
            L = tpmm2t[j][1];

            e = TPM::gt2s(K,0);
            z = TPM::gt2s(K,1);

            t = TPM::gt2s(L,0);
            h = TPM::gt2s(L,1);

            (*this)(i,j) = 0.0;

            if(z <= t){

               for(int k = 0;k < a;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(k,a,b,t,l,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(k,a,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = a + 1;k < c;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = c + 1;k < b;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,k,d,l,e,z);

                  for(int l = e + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = z + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,k,d,e,z,l);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,k,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(c,k,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,k,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(c,k,d,e,z,l);

               }

               for(int k = b + 1;k < d;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,k,d,l,e,z);

                  for(int l = e + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = z + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,k,d,e,z,l);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,k,d,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,k,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,k,d,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,k,d,e,z,l);

               }

               for(int k = d + 1;k < M;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,l,e,z);

                  for(int l = e + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = z + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,e,z,l);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,d,k,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,d,k,e,z,l);

               }

            }
            else if(z <= h){

               for(int k = 0;k < a;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = t + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(k,a,b,t,l,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(k,a,b,t,l,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(k,a,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = a + 1;k < c;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = t + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = c + 1;k < b;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,k,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = t + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,k,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = z + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,k,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(c,k,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,k,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(c,k,d,e,z,l);

               }

               for(int k = b + 1;k < d;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,k,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = t + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,k,d,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = z + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,k,d,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,k,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,k,d,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,k,d,e,z,l);

               }

               for(int k = d + 1;k < M;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = t + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = z + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,d,k,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,d,k,e,z,l);

               }

            }
            else{

               for(int k = 0;k < a;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(k,a,b,t,l,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = h + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(k,a,b,t,h,l) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(k,a,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = a + 1;k < c;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = h + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = c + 1;k < b;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,k,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,k,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = h + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,k,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = z + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,k,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(c,k,d,e,z,l);

               }

               for(int k = b + 1;k < d;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,k,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,k,d,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = h + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,k,d,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = z + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,k,d,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,k,d,e,z,l);

               }

               for(int k = d + 1;k < M;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = h + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = z + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,d,k,e,z,l);

               }

            }
         
         }

      }
      else{

         for(int j = i;j < gn();++j){

            K = tpmm2t[j][0];
            L = tpmm2t[j][1];

            e = TPM::gt2s(K,0);
            z = TPM::gt2s(K,1);

            t = TPM::gt2s(L,0);
            h = TPM::gt2s(L,1);

            (*this)(i,j) = 0.0;

            if(z <= t){

               for(int k = 0;k < a;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(k,a,b,t,l,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(k,a,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = a + 1;k < c;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = c + 1;k < d;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,k,d,l,e,z);

                  for(int l = e + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = z + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,k,d,e,z,l);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,k,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(c,k,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,k,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(c,k,d,e,z,l);

               }

               for(int k = d + 1;k < b;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,d,k,l,e,z);

                  for(int l = e + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = z + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,d,k,e,z,l);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,d,k,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(c,d,k,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,d,k,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(c,d,k,e,z,l);

               }

               for(int k = b + 1;k < M;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,l,e,z);

                  for(int l = e + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = z + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,e,z,l);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,d,k,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,d,k,e,z,l);

               }

            }
            else if(z <= h){

               for(int k = 0;k < a;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = t + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(k,a,b,t,l,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(k,a,b,t,l,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(k,a,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = a + 1;k < c;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = t + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(k,c,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = c + 1;k < d;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,k,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = t + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,k,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = z + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,k,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(c,k,d,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,k,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(c,k,d,e,z,l);

               }

               for(int k = d + 1;k < b;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,d,k,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = t + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,d,k,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = z + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,d,k,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(c,d,k,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,d,k,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(c,d,k,e,z,l);

               }

               for(int k = b + 1;k < M;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = t + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = z + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,d,k,e,z,l);

                  for(int l = h + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,d,k,e,z,l);

               }

            }
            else{

               for(int k = 0;k < a;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(k,a,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(k,a,b,t,l,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = h + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(k,a,b,e,l,z) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(k,a,b,t,h,l) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(k,a,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(k,a,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = a + 1;k < c;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = h + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(k,c,d,e,l,z);

                  for(int l = z + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(k,c,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(k,c,d,e,z,l);

               }

               for(int k = c + 1;k < d;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,k,d,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,k,d,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,k,d,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = h + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,k,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(c,k,d,e,l,z);

                  for(int l = z + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,k,d,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(c,k,d,e,z,l);

               }

               for(int k = d + 1;k < b;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,l,e,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,d,k,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,k,b,l,t,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,d,k,t,l,h) + dpm.ordacc(a,k,b,t,l,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = h + 1;l < z;++l)
                     (*this)(i,j) += dpm.ordacc(a,k,b,e,l,z) * dpm.ordacc(c,d,k,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = z + 1;l < M;++l)
                     (*this)(i,j) -= dpm.ordacc(a,k,b,e,z,l) * dpm.ordacc(c,d,k,t,h,l) + dpm.ordacc(a,k,b,t,h,l) * dpm.ordacc(c,d,k,e,z,l);

               }

               for(int k = b + 1;k < M;++k){

                  for(int l = 0;l < e;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,l,e,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,l,e,z);

                  for(int l = e + 1;l < t;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,l,t,h) + dpm.ordacc(a,b,k,l,t,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = t + 1;l < h;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,t,l,h) + dpm.ordacc(a,b,k,t,l,h) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = h + 1;l < z;++l)
                     (*this)(i,j) -= dpm.ordacc(a,b,k,e,l,z) * dpm.ordacc(c,d,k,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,d,k,e,l,z);

                  for(int l = z + 1;l < M;++l)
                     (*this)(i,j) += dpm.ordacc(a,b,k,e,z,l) * dpm.ordacc(c,d,k,t,h,l) + dpm.ordacc(a,b,k,t,h,l) * dpm.ordacc(c,d,k,e,z,l);

               }

            }
         
         }

      }

   }

   this->symmetrize();

}
