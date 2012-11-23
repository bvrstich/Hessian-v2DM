#ifndef HESSIAN_H
#define HESSIAN_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

/**
 * @author Brecht Verstichel
 * @date 23-11-2012\n\n
 * This class Hessian is a class written for matrices of two particle matrices, it inherits alle the function from its mother 
 * Matrix, some special member functions and two lists that give the relationship between the sp and the tp 
 * basis.
 */

class Hessian : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param hess_p the Hessian you want to print
    */
   friend ostream &operator<<(ostream &output,const Hessian &hess_p);

   public:
      
      //constructor
      Hessian();

      //copy constructor
      Hessian(const Hessian &);

      //destructor
      virtual ~Hessian();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d,int e,int z,int t,int h) const;

      //access to the numbers in tp mode
      double operator()(int I,int J,int K,int L) const;

      void I(const TPM &);

      void Q(const TPM &);

      void lagr();

      static int gn();

      static int ghess2t(int,int);

      static int gt2hess(int,int);
      
      static void init();

      static void clear();

   private:

      //!list relating the single-particle space to the Hessian basis
      static vector< vector<int> > hess2t;

      //!list relating the single-particle space to the Hessian basis
      static int **t2hess;

};

#endif
