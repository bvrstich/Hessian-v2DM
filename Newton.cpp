#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

vector< vector<int> > Newton::hess2t;
int **Newton::t2hess;

vector<double> Newton::norm;

/**
 * allocate and fill the static variables and lists
 */
void Newton::init(){

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

         if(i == j)
            norm.push_back(std::sqrt(0.5));
         else
            norm.push_back(1.0);

         ++hess;

      }

}

/**
 * clear the static lists
 */
void Newton::clear(){

   for(int i = 0;i < TPM::gn();++i)
      delete [] t2hess[i];

   delete [] t2hess;

}

/**
 * standard constructor
 */
Newton::Newton(){

   H = new Matrix(hess2t.size());

   x = new double [hess2t.size()];

}

/**
 * copy constructor
 */
Newton::Newton(const Newton &newton_c){

   H = new Matrix(newton_c.gH());

   x = new double [hess2t.size()];

   for(unsigned int i = 0;i < hess2t.size();++i)
      x[i] = newton_c.gx()[i];

}

/**
 * destructor
 */
Newton::~Newton(){

   delete H;

   delete [] x;

}

/**
 * @return dimension of the hessian matrix
 */
int Newton::gdim() const {

   return hess2t.size();

}

/**
 * access to the lists from outside the class
 */
int Newton::gt2hess(int i,int j){

   return t2hess[i][j];

}

/**
 * access to the lists from outside the class
 * @param option == 0 return a, == 1 return b
 */
int Newton::ghess2t(int i,int option){

   return hess2t[i][option];

}

/**
 * @return the full Hessian matrix, read only
 */
const Matrix &Newton::gH() const {

   return *H;

}

/**
 * @return the "vector" x, read only
 */
const double* Newton::gx() const {

   return x;

}

/**
 * @return the "vector" x, read and write
 */
double* Newton::gx() {

   return x;

}

/**
 * construct the different parts of the linear system
 * @param t potential scaling parameter
 * @param ham hamiltonian object
 * @param P object containing the inverse of the matrix constraints p and q.
 */
void Newton::construct(double t,const TPM &ham,const SUP &P){

   //first construct the gradient
   constr_grad(t,ham,P);

   //construct the p part of the hessian
   constr_hess_I(t,P);

   H->symmetrize();

   //invert the hessian using cholesky
   H->invert();

   //construct the matrix-vector product of inverse hessian with gradient
   double *y = new double [hess2t.size()];

   H->mv(1.0,x,y);

   //now we can obtain the lagrange multiplier
   double lambda = vectrace(y)/Htrace();

   //construct the half trace inverse hessian:
   Hbar(lambda,x);

   //now we can construct the solution:
   for(unsigned int i = 0;i < hess2t.size();++i)
      x[i] = y[i] - x[i];

   delete [] y;

}

/**
 * construct 'minus' the gradient vector
 * @param t potential scaling parameter
 * @param ham hamiltonian object
 * @param P object containing the inverse of the matrix constraints p and q.
 */
void Newton::constr_grad(double t,const TPM &ham,const SUP &P){

   int I,J;

   for(unsigned int i = 0;i < hess2t.size();++i){

      I = hess2t[i][0];
      J = hess2t[i][1];

      x[i] = t * P.gI()(I,J) - ham(I,J);

      x[i] *= std::sqrt(2.0) * norm[i];

   }

}

/**
 * access to the norm from outside of the class
 * @param i hess index, if a == b norm = 1.0/sqrt(2.0)
 */
double Newton::gnorm(int i){

   return norm[i];

}

/**
 * construct the I part of the hessian matrix
 */
void Newton::constr_hess_I(double t,const SUP &P){

   int I,J,K,L;

   for(unsigned int i = 0;i < hess2t.size();++i){

      I = hess2t[i][0];
      J = hess2t[i][1];

      for(unsigned int j = i;j < hess2t.size();++j){

         K = hess2t[j][0];
         L = hess2t[j][1];

         (*H)(i,j) = ( P.gI()(I,K) * P.gI()(J,L) + P.gI()(I,L) * P.gI()(J,K) ) * norm[i] * norm[j] * t;

      }

   }

}

/**
 * access to the x-elements, the solution, as it were.
 * elements are already transformed to matrix mode
 */
double Newton::gx(int I,int J) const {

   int i = t2hess[I][J];

   return x[i] / (std::sqrt(2.0) * norm[i]);

}

/**
 * trace a vector as if it were a matrix, using the list t2hess
 * @param y input vector
 * @return the trace
 */
double Newton::vectrace(const double *y) const {

   double ward = 0.0;

   for(int I = 0;I < TPM::gn();++I)
      ward += y[t2hess[I][I]];

   return ward;

}

/**
 * @return the trace of the hessian, with trace I actually mean \sum_{IJ} H_{II;JJ}
 */
double Newton::Htrace() const {

   double ward = 0.0;

   for(int I = 0;I < TPM::gn();++I)
      for(int J = 0;J < TPM::gn();++J)
      ward += (*H)(t2hess[I][I],t2hess[J][J]);

   return ward;



}

/**
 * map the hessian onto a vector by tracing out the rows or columns
 */
void Newton::Hbar(double scale,double *y){

   for(unsigned int i = 0;i < hess2t.size();++i){

      y[i] = 0.0;

      for(int K = 0;K < TPM::gn();++K)
         y[i] += (*H)(i,t2hess[K][K]);

      y[i] *= scale;

   }

}
