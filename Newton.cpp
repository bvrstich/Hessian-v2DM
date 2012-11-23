#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

double *Newton::norm;

/**
 * allocate and fill the static variables and lists
 */
void Newton::init(){

   norm = new double [Hessian::gn()];

   int hess = 0;

   for(int i = 0;i < TPM::gn();++i)
      for(int j = i;j < TPM::gn();++j){

         if(i == j)
            norm[hess] = 0.5;
         else
            norm[hess] = std::sqrt(0.5);

         ++hess;

      }

}

/**
 * clear the static lists
 */
void Newton::clear(){

   delete [] norm;

}

/**
 * standard constructor
 */
Newton::Newton(){

   H = new Hessian();

   x = new double [Hessian::gn() + 1];

}

/**
 * copy constructor
 */
Newton::Newton(const Newton &newton_c){

   H = new Hessian(newton_c.gH());

   x = new double [Hessian::gn() + 1];

   for(int i = 0;i < Hessian::gn() + 1;++i)
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
 * @return the full Hessian matrix, read only
 */
const Hessian &Newton::gH() const {

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
   gradient(t,ham,P);

   //construct the p part of the hessian
   H->I(P.gI());

#ifdef __Q_CON
   H->Q(P.gQ());
#endif

   H->dscal(t);

   //the constraint/lagrange multiplier part of the Hessian
   H->lagr();

   //last element zero!
   (*H)(Hessian::gn(),Hessian::gn()) = 0.0;
   
   H->symmetrize();

   //and last but not least, solve the system
   H->solve_sy(x);

}

/**
 * construct 'minus' the gradient vector
 * @param t potential scaling parameter
 * @param ham hamiltonian object
 * @param P object containing the inverse of the matrix constraints p and q.
 */
void Newton::gradient(double t,const TPM &ham,const SUP &P){

   int I,J;

#ifdef __Q_CON
   TPM QQ;
   QQ.Q(1,P.gQ());
#endif

   for(int i = 0;i < Hessian::gn();++i){

      I = Hessian::ghess2t(i,0);
      J = Hessian::ghess2t(i,1);

      x[i] = t * P.gI()(I,J) - ham(I,J);

#ifdef __Q_CON
      x[i] += t * QQ(I,J);
#endif

      x[i] *= 2.0 * norm[i];

   }

   //last part of right-hand side (lagrange multiplier)
   x[Hessian::gn()] = 0.0;

}

/**
 * access to the norm from outside of the class
 * @param i hess index, if a == b norm = 1.0/sqrt(2.0)
 */
double Newton::gnorm(int i){

   return norm[i];

}

/**
 * access to the norm from outside of the class, in tp mode
 * @param I row index
 * @param J column index
 */
double Newton::gnorm(int I,int J){

   return norm[Hessian::gt2hess(I,J)];

}

/**
 * access to the x-elements, the solution, as it were.
 * elements are already transformed to matrix mode
 */
double Newton::gx(int I,int J) const {

   int i = Hessian::gt2hess(I,J);

   return x[i];

}

/**
 * access to the x-elements, the solution, as it were in sp mode
 * elements are already transformed to matrix mode
 */
double Newton::gx(int a,int b,int c,int d) const {

   if( (a == b) || (c == d) )
      return 0.0;

   int phase = 1;

   if(a > b)
      phase *= -1;

   if(c > d)
      phase *= -1;

   int I = TPM::gs2t(a,b);
   int J = TPM::gs2t(c,d);

   int i = Hessian::gt2hess(I,J);

   return phase * x[i];

}
