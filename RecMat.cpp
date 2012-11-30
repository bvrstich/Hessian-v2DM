#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "RecMat.h"
#include "lapack.h"

/**
 * constructor 
 * @param m number of rows
 * @param n number of columns
 */
RecMat::RecMat(int m,int n){

   this->m = m;
   this->n = n;

   recmat = new double * [n];
   recmat[0] = new double [n*m];

   for(int i = 1;i < n;++i)
      recmat[i] = recmat[i - 1] + m;

}

/**
 * copy constructor 
 * @param rm_copy The recmat you want to be copied into the object you are constructing
 */
RecMat::RecMat(const RecMat &rm_copy){

   this->m = rm_copy.gm();
   this->n = rm_copy.gn();

   recmat = new double * [n];
   recmat[0] = new double [n*m];

   for(int i = 1;i < n;++i)
      recmat[i] = recmat[i - 1] + m;

   int dim = n*m;
   int incx = 1;
   int incy = 1;

   dcopy_(&dim,rm_copy.recmat[0],&incx,recmat[0],&incy);

}

/**
 * Destructor
 */
RecMat::~RecMat(){

   delete [] recmat[0];
   delete [] recmat;

}

/**
 * overload the equality operator
 * @param recrm_copy The recmat you want to be copied into this
 */
RecMat &RecMat::operator=(const RecMat &recrm_copy){

   int dim = n*m;
   int incx = 1;
   int incy = 1;

   dcopy_(&dim,recrm_copy.recmat[0],&incx,recmat[0],&incy);

   return *this;

}

/**
 * Make all the number in your recmat equal to the number a, e.g. usefull for initialization (RecMat M = 0)
 * @param a the number
 */
RecMat &RecMat::operator=(double a){

   for(int i = 0;i < m;++i)
      for(int j = 0;j < n;++j)
         recmat[j][i] = a;

   return *this;

}

/**
 * overload the += operator for matrices
 * @param recmat_pl The recmat you want to add to this
 */
RecMat &RecMat::operator+=(const RecMat &recmat_pl){

   int dim = n*m;
   int inc = 1;
   double alpha = 1.0;

   daxpy_(&dim,&alpha,recmat_pl.recmat[0],&inc,recmat[0],&inc);

   return *this;

}

/**
 * overload the -= operator for matrices
 * @param recmat_pl The recmat you want to deduct from this
 */
RecMat &RecMat::operator-=(const RecMat &recmat_pl){

   int dim = n*m;
   int inc = 1;
   double alpha = -1.0;

   daxpy_(&dim,&alpha,recmat_pl.recmat[0],&inc,recmat[0],&inc);

   return *this;

}

/**
 * add the recmat recmat_pl times the constant alpha to this
 * @param alpha the constant to multiply the recmat_pl with
 * @param recmat_pl the RecMat to be multiplied by alpha and added to this
 */
RecMat &RecMat::daxpy(double alpha,const RecMat &recmat_pl){

   int dim = n*m;
   int inc = 1;

   daxpy_(&dim,&alpha,recmat_pl.recmat[0],&inc,recmat[0],&inc);

   return *this;

}
/**
 * /= operator overloaded: divide by a constant
 * @param c the number to divide your recmat through
 */
RecMat &RecMat::operator/=(double c){

   int dim = n*m;
   int inc = 1;

   double alpha = 1.0/c;

   dscal_(&dim,&alpha,recmat[0],&inc);

   return *this;

}

/**
 * write access to your recmat, change the number on row i and column j
 * remark that for the conversion to lapack functions the double pointer is transposed!
 * @param i row number
 * @param j column number
 * @return the entry on place i,j
 */
double &RecMat::operator()(int i,int j){

   return recmat[j][i];

}

/**
 * read access to your recmat, view the number on row i and column j
 * remark that for the conversion to lapack functions the double pointer is transposed!
 * @param i row number
 * @param j column number
 * @return the entry on place i,j
 */
double RecMat::operator()(int i,int j) const {

   return recmat[j][i];

}

/**
 * @return the underlying pointer to recmat, useful for mkl applications
 */
double **RecMat::gRecMat(){

   return recmat;

}

/**
 * const version
 * @return the underlying pointer to recmat, useful for mkl applications
 */
double **RecMat::gRecMat() const {

   return recmat;

}

/**
 * @return the number of rows
 */
int RecMat::gm() const{

   return m;

}

/**
 * @return the number of columns
 */
int RecMat::gn() const{

   return n;

}

/**
 * Scale the recmat (*this) with parameter alpha
 * @param alpha scalefactor
 */
void RecMat::dscal(double alpha){

   int dim = n*m;
   int inc = 1;

   dscal_(&dim,&alpha,recmat[0],&inc);

}

/**
 * Fill the recmat with random numbers.
 */
void RecMat::fill_Random(){

   srand(time(NULL));

   for(int i = 0;i < m;++i)
      for(int j = 0;j < n;++j)
         recmat[j][i] = (double) rand()/RAND_MAX;

}

ostream &operator<<(ostream &output,const RecMat &recmat_p){

   for(int i = 0;i < recmat_p.gm();++i)
      for(int j = 0;j < recmat_p.gn();++j)
         output << i << "\t" << j << "\t" << recmat_p(i,j) << endl;

   return output;

}
