#ifndef RECMAT_H
#define RECMAT_H

#include <iostream>
#include <cstdlib>

using std::ostream;

#include "Vector.h"

/**
 * @author Brecht Verstichel
 * @date 18-02-2010\n\n
 * This is a class written for symmetric matrices. It is a wrapper around a double pointer and
 * redefines much used lapack and blas routines as memberfunctions
 */

class RecMat{

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param rm_p de RecMat you want to print
    */
   friend ostream &operator<<(ostream &output,const RecMat &rm_p);

   public:

      //constructor
      RecMat(int m,int n);

      //copy constructor
      RecMat(const RecMat &);

      //construct with filename
      RecMat(const char *filename);

      //destructor
      virtual ~RecMat();

      //overload equality operator
      RecMat &operator=(const RecMat &);

      RecMat &operator=(double );

      //overload += operator
      RecMat &operator+=(const RecMat &);

      //overload -= operator
      RecMat &operator-=(const RecMat &);

      RecMat &daxpy(double alpha,const RecMat &);

      RecMat &operator/=(double );

      RecMat &mprod(const RecMat &,const RecMat &);

      //easy to change the numbers
      double &operator()(int i,int j);

      //easy to access the numbers
      double operator()(int i,int j) const;

      //get the pointer to the matrix
      double **gRecMat();

      //get the pointer to the matrix
      double **gRecMat() const;

      int gn() const;
      
      int gm() const;

      void dscal(double alpha);

      void fill_Random();

   private:

      //!double pointer of doubles, contains the numbers, the matrix
      double **recmat;

      //!number of rows
      int m;

      //!number of columns
      int n;

};

#endif
