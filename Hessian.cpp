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
Hessian::Hessian() : Matrix(TPTPM::gn() + 1) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix hess_c
 * @param hess_c object that will be copied into this.
 */
Hessian::Hessian(const Hessian &hess_c) : Matrix(hess_c){ }

/**
 * destructor
 */
Hessian::~Hessian(){ }

ostream &operator<<(ostream &output,const Hessian &hess_p){

   int I,J,K,L;

   int a,b,c,d;
   int e,z,t,h;

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

         output << i << "\t" << j << "\t|\t" << I << "\t" << J << "\t" << K << "\t" << L << "\t|\t" << 
         
            "(" << a << "," << b << "," << c << "," << d << ")\t(" << e << "," << z << "," << t << "," << h << ")\t|\t" << hess_p(i,j) << endl;

      }

   }

   return output;

}

/**
 * construct the I part of the hessian matrix
 */
void Hessian::I(const TPM &tpm){

   int I,J,K,L;

   for(int i = 0;i < TPTPM::gn();++i){

      I = TPTPM::gtpmm2t(i,0);
      J = TPTPM::gtpmm2t(i,1);

      for(int j = i;j < TPTPM::gn();++j){

         K = TPTPM::gtpmm2t(j,0);
         L = TPTPM::gtpmm2t(j,1);

         (*this)(i,j) = 2.0 * ( tpm(I,K) * tpm(J,L) + tpm(I,L) * tpm(J,K) ) * Newton::gnorm(i) * Newton::gnorm(j);

      }

   }

}

/**
 * construct the lagrange multiplier part of the Hessian
 */
void Hessian::lagr(){

   int I,J;

   for(int i = 0;i < TPTPM::gn();++i){

      I = TPTPM::gtpmm2t(i,0);
      J = TPTPM::gtpmm2t(i,1);

      if(I == J)
         (*this)(i,TPTPM::gn()) = 1.0;
      else
         (*this)(i,TPTPM::gn()) = 0.0;

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

/**
 * construct the G part of the hessian matrix: add to current hessian!
 */
void Hessian::G(const PHM &G){

   PHSPM phspm;
   phspm.dirprodtrace(2.0/(Tools::gN() - 1.0),G);

   SPSPM spspm;
   spspm.bar(1.0/(Tools::gN() - 1.0),phspm);

   int I,J,K,L;

   int a,b,c,d;
   int e,z,t,h;

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
         (*this)(i,j) += 2.0 * Newton::gnorm(i) * Newton::gnorm(j) * (

               G(a,d,e,h) * G(c,b,t,z) + G(a,d,t,z) * G(c,b,e,h) - G(a,d,z,h) * G(c,b,t,e) - G(a,d,t,e) * G(c,b,z,h)

               - G(a,d,e,t) * G(c,b,h,z) - G(a,d,h,z) * G(c,b,e,t) + G(a,d,z,t) * G(c,b,h,e) + G(a,d,h,e) * G(c,b,z,t)

               - G(b,d,e,h) * G(c,a,t,z) - G(b,d,t,z) * G(c,a,e,h) + G(b,d,z,h) * G(c,a,t,e) + G(b,d,t,e) * G(c,a,z,h)

               + G(b,d,e,t) * G(c,a,h,z) + G(b,d,h,z) * G(c,a,e,t) - G(b,d,z,t) * G(c,a,h,e) - G(b,d,h,e) * G(c,a,z,t)

               - G(a,c,e,h) * G(d,b,t,z) - G(a,c,t,z) * G(d,b,e,h) + G(a,c,z,h) * G(d,b,t,e) + G(a,c,t,e) * G(d,b,z,h)

               + G(a,c,e,t) * G(d,b,h,z) + G(a,c,h,z) * G(d,b,e,t) - G(a,c,z,t) * G(d,b,h,e) - G(a,c,h,e) * G(d,b,z,t)

               + G(b,c,e,h) * G(d,a,t,z) + G(b,c,t,z) * G(d,a,e,h) - G(b,c,z,h) * G(d,a,t,e) - G(b,c,t,e) * G(d,a,z,h)

               - G(b,c,e,t) * G(d,a,h,z) - G(b,c,h,z) * G(d,a,e,t) + G(b,c,z,t) * G(d,a,h,e) + G(b,c,h,e) * G(d,a,z,t) );

         if(b == d){

            //first four PHSPM terms:
            (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * ( phspm(e,h,t,z,a,c) - phspm(z,h,t,e,a,c) - phspm(e,t,h,z,a,c) + phspm(z,t,h,e,a,c) );

            //then four SPSPM terms:
            if(z == h)
               (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * spspm(a,c,e,t);

            if(z == t)
               (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * spspm(a,c,e,h);

            if(e == t)
               (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * spspm(a,c,z,h);

         }

         if(b == c){

            //first four PHSPM terms:
            (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * ( phspm(e,h,t,z,a,d) - phspm(z,h,t,e,a,d) - phspm(e,t,h,z,a,d) + phspm(z,t,h,e,a,d) );

            //then four SPSPM terms
            if(z == h)
               (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * spspm(a,d,e,t);

            if(z == t)
               (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * spspm(a,d,e,h);

            if(e == t)
               (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * spspm(a,d,z,h);

         }

         if(a == c){

            //first four PHSPM terms:
            (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * ( phspm(e,h,t,z,b,d) - phspm(z,h,t,e,b,d) - phspm(e,t,h,z,b,d) + phspm(z,t,h,e,b,d) );

            //then four SPSPM terms
            if(z == h)
               (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * spspm(b,d,e,t);

            if(z == t)
               (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * spspm(b,d,e,h);

            if(e == t)
               (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * spspm(b,d,z,h);

         }

         if(z == h)
            (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * ( phspm(a,d,c,b,e,t) - phspm(b,d,c,a,e,t) - phspm(a,c,d,b,e,t) + phspm(b,c,d,a,e,t) );

         if(z == t)
            (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * ( phspm(a,d,c,b,e,h) - phspm(b,d,c,a,e,h) - phspm(a,c,d,b,e,h) + phspm(b,c,d,a,e,h) );

         if(e == t)
            (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * ( phspm(a,d,c,b,z,h) - phspm(b,d,c,a,z,h) - phspm(a,c,d,b,z,h) + phspm(b,c,d,a,z,h) );

      }

   }

}

/**
 * construct the T1 part of the hessian matrix: add to current hessian!
 * @param dpm input DPM object
 */
void Hessian::T(const DPM &dpm){

   TPTPM tpmm;
   tpmm.dirprodtrace(dpm);

   TPSPM tpspm;
   tpspm.bar(1.0/(Tools::gN() - 1.0),tpmm);

   SPSPM spmm;
   spmm.bar(0.5/(Tools::gN() - 1.0),tpspm);

   DPM T2;
   T2.squaresym(dpm);

   TPM T2bar;
   T2bar.bar(4.0/(Tools::gN() * (Tools::gN() - 1.0)),T2);

   SPM T2barbar;
   T2barbar.bar(0.5/(Tools::gN() - 1.0),T2bar);

   double trace = 4.0/ ( Tools::gN() * Tools::gN() * (Tools::gN() - 1.0) * (Tools::gN() - 1.0) ) * T2.trace();

   int I,J,K,L;

   int a,b,c,d;
   int e,z,t,h;

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

         //first direct product term:
         (*this)(i,j) += 2.0 * Newton::gnorm(i) * Newton::gnorm(j) * tpmm(i,j);

         if(I == J){

            //np
            if(K == L)
               (*this)(i,j) += trace;

            //tp
            (*this)(i,j) += Newton::gnorm(j) * T2bar(e,z,t,h);

            //4 sp
            if(z == h)
               (*this)(i,j) -= Newton::gnorm(j) * T2barbar(e,t);

            if(z == t)
               (*this)(i,j) += Newton::gnorm(j) * T2barbar(e,h);

            if(e == t)
               (*this)(i,j) -= Newton::gnorm(j) * T2barbar(z,h);

         }

         if(K == L){

            //tp
            (*this)(i,j) += Newton::gnorm(i) * T2bar(a,b,c,d);

            //4 sp
            if(b == d)
               (*this)(i,j) -= Newton::gnorm(i) * T2barbar(a,c);

            if(b == c)
               (*this)(i,j) += Newton::gnorm(i) * T2barbar(a,d);

            if(a == c)
               (*this)(i,j) -= Newton::gnorm(i) * T2barbar(b,d);

         }

         if(b == d){

            (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * tpspm(e,z,t,h,a,c);

            if(z == h)
               (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * spmm(a,c,e,t);

            if(z == t)
               (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * spmm(a,c,e,h);

            if(e == t)
               (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * spmm(a,c,z,h);

         }

         if(b == c){

            (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * tpspm(e,z,t,h,a,d);

            if(z == h)
               (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * spmm(a,d,e,t);

            if(z == t)
               (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * spmm(a,d,e,h);

            if(e == t)
               (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * spmm(a,d,z,h);

         }

         if(a == c){

            (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * tpspm(e,z,t,h,b,d);

            if(z == h)
               (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * spmm(b,d,e,t);

            if(z == t)
               (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * spmm(b,d,e,h);

            if(e == t)
               (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * spmm(b,d,z,h);

         }

         if(z == h)
            (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * tpspm(a,b,c,d,e,t);

         if(z == t)
            (*this)(i,j) += Newton::gnorm(i) * Newton::gnorm(j) * tpspm(a,b,c,d,e,h);

         if(e == t)
            (*this)(i,j) -= Newton::gnorm(i) * Newton::gnorm(j) * tpspm(a,b,c,d,z,h);

      }
   }

   this->symmetrize();

}
