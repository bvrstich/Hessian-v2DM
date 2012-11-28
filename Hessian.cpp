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

int **Hessian::t2hess;
vector< vector<int> > Hessian::hess2t;

/**
 * initialize the static lists
 */
void Hessian::init(){

   t2hess = new int * [TPM::gn()];

   for(int i = 0;i < TPM::gn();++i)
      t2hess[i] = new int [TPM::gn()];

   vector<int> v(2);

   int hess = 0;

   for(int i = 0;i < TPM::gn();++i)
      for(int j = i;j < TPM::gn();++j){

         v[0] = i;
         v[1] = j;

         hess2t.push_back(v);

         t2hess[i][j] = hess;
         t2hess[j][i] = hess;

         ++hess;

      }

}

/**
 * deallocate the static lists
 */
void Hessian::clear(){

   for(int i = 0;i < TPM::gn();++i)
      delete [] t2hess[i];

   delete [] t2hess;

}

/**
 * standard constructor:
 */
Hessian::Hessian() : Matrix(hess2t.size() + 1) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix hess_c
 * @param hess_c object that will be copied into this.
 */
Hessian::Hessian(const Hessian &hess_c) : Matrix(hess_c){ }

/**
 * destructor
 */
Hessian::~Hessian(){ }

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
double Hessian::operator()(int a,int b,int c,int d,int e,int z,int t,int h) const{

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

      int i = t2hess[I][J];
      int j = t2hess[K][L];

      return phase * (*this)(i,j);

   }

}

/**
 * access the elements of the matrix in tp mode
 * @param I first tp index that forms the hess row index i together with J
 * @param J second tp index that forms the hess row index i together with I
 * @param K first tp index that forms the hess column index j together with L
 * @param L second tp index that forms the hess column index j together with K
 * @return the number on place Hessian(i,j)
 */
double Hessian::operator()(int I,int J,int K,int L) const{

   int i = t2hess[I][J];
   int j = t2hess[K][L];

   return (*this)(i,j);

}

ostream &operator<<(ostream &output,const Hessian &hess_p){

   for(int i = 0;i < hess_p.gn();++i)
      for(int j = 0;j < hess_p.gn();++j){

         output << i << "\t" << j << "\t|\t" << hess_p.hess2t[i][0] << "\t" << hess_p.hess2t[i][1]

            << "\t" << hess_p.hess2t[j][0] << "\t" << hess_p.hess2t[j][1] << "\t" << hess_p(i,j) << endl;

      }

   return output;

}

/**
 * @return the dimension of a Hessian matrix
 */
int Hessian::gn(){

   return hess2t.size();

}

/**
 * access to the lists from outside the class
 */
int Hessian::gt2hess(int i,int j){

   return t2hess[i][j];

}

/**
 * access to the lists from outside the class
 * @param option == 0 return a, == 1 return b
 */
int Hessian::ghess2t(int i,int option){

   return hess2t[i][option];

}

/**
 * construct the I part of the hessian matrix
 */
void Hessian::I(const TPM &tpm){

   int I,J,K,L;

   for(unsigned int i = 0;i < hess2t.size();++i){

      I = hess2t[i][0];
      J = hess2t[i][1];

      for(unsigned int j = i;j < hess2t.size();++j){

         K = hess2t[j][0];
         L = hess2t[j][1];

         (*this)(i,j) = 2.0 * ( tpm(I,K) * tpm(J,L) + tpm(I,L) * tpm(J,K) ) * Newton::gnorm(i) * Newton::gnorm(j);

      }

   }

}

/**
 * construct the lagrange multiplier part of the Hessian
 */
void Hessian::lagr(){

   int I,J;

   for(unsigned int i = 0;i < hess2t.size();++i){

      I = hess2t[i][0];
      J = hess2t[i][1];

      if(I == J)
         (*this)(i,hess2t.size()) = 1.0;
      else
         (*this)(i,hess2t.size()) = 0.0;

   }

}

/**
 * construct the Q part of the hessian matrix: add to current hessian!
 */
void Hessian::Q(const TPM &Q){

   TPSPM tpspm;
   tpspm.dirprodtrace(2.0/(Tools::gN() - 1.0),Q);

   SPSPM spspm;
   spspm.bar(1.0/(Tools::gN() - 1.0),tpspm);

   TPM Q2;
   Q2.squaresym(Q);

   SPM Q2bar;
   Q2bar.bar(4.0/(Tools::gN()*(Tools::gN()-1.0)*(Tools::gN()-1.0)),Q2);

   double trace = Q2bar.trace()/(double)Tools::gN();

   int I,J,K,L;

   int a,b,c,d;
   int e,z,t,h;

   for(unsigned int i = 0;i < hess2t.size();++i){

      I = hess2t[i][0];
      J = hess2t[i][1];

      a = TPM::gt2s(I,0);
      b = TPM::gt2s(I,1);

      c = TPM::gt2s(J,0);
      d = TPM::gt2s(J,1);

      for(unsigned int j = i;j < hess2t.size();++j){

         K = hess2t[j][0];
         L = hess2t[j][1];

         e = TPM::gt2s(K,0);
         z = TPM::gt2s(K,1);

         t = TPM::gt2s(L,0);
         h = TPM::gt2s(L,1);

         //first direct product term
         (*this)(i,j) += 2.0 * Newton::gnorm(i) * Newton::gnorm(j) * ( Q(I,K) * Q(J,L) + Q(I,L) * Q(J,K) );

         //first the diagonal terms on tp space
         if(I == J){

            if(K == L)
               (*this)(i,j) += trace;

            (*this)(i,j) += 4.0 / ( Tools::gN() * (Tools::gN() - 1.0) ) * Newton::gnorm(j) * Q2(K,L);

            if(e == t)
               (*this)(i,j) -= Newton::gnorm(j) * Q2bar(z,h);
            
            if(z == t)
               (*this)(i,j) += Newton::gnorm(j) * Q2bar(e,h);

            if(z == h)
               (*this)(i,j) -= Newton::gnorm(j) * Q2bar(e,t);

         }

         if(K == L){

            (*this)(i,j) += 4.0 / ( Tools::gN() * (Tools::gN() - 1.0) ) * Newton::gnorm(i) * Q2(I,J);

            if(a == c)
               (*this)(i,j) -= Newton::gnorm(i) * Q2bar(b,d);

            if(b == c)
               (*this)(i,j) += Newton::gnorm(i) * Q2bar(a,d);

            if(b == d)
               (*this)(i,j) -= Newton::gnorm(i) * Q2bar(a,c);

         }

         //Now the equalities in single-particle indices
         if(a == c){

            (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * tpspm(e,z,t,h,b,d);

            //4 spspm terms
            if(e == t)
               (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * spspm(b,d,z,h);

            if(z == t)
               (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * spspm(b,d,e,h);

            if(z == h)
               (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * spspm(b,d,e,t);

         }

         if(b == c){

            (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * tpspm(e,z,t,h,a,d);

            //4 spspm terms
            if(e == t)
               (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * spspm(a,d,z,h);

            if(z == t)
               (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * spspm(a,d,e,h);

            if(z == h)
               (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * spspm(a,d,e,t);

         }

         if(b == d){

            (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * tpspm(e,z,t,h,a,c);

            //4 spspm terms
            if(e == t)
               (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * spspm(a,c,z,h);

            if(z == t)
               (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * spspm(a,c,e,h);

            if(z == h)
               (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * spspm(a,c,e,t);

         }

         //four more years! I mean terms
         if(e == t)
            (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * tpspm(a,b,c,d,z,h);

         if(z == t)
            (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * tpspm(a,b,c,d,e,h);

         if(z == h)
            (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * tpspm(a,b,c,d,e,t);

      }

   }

}
