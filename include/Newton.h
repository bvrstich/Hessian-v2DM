#ifndef NEWTON_H
#define NEWTON_H

#include <iostream>
#include <cstdlib>
#include <vector>

using std::ostream;
using std::vector;

#include "Vector.h"

/**
 * @author Brecht Verstichel
 * @date 21-11-2012\n\n
 * This is a class written for the construction and solution of the Newton system. 
 */

class Newton{

   public:

      //constructor
      Newton();

      //copy constructor
      Newton(const Newton &);

      //construct with filename
      Newton(const char *filename);

      //destructor
      virtual ~Newton();

      int gdim() const;

      const Matrix &gH() const;

      const double *gx() const;

      double *gx();

      double gx(int,int) const;

      double vectrace(const double *) const;

      double Htrace() const;

      void Hbar(double,double *);

      void construct(double,const TPM &,const SUP &);

      void constr_grad(double,const TPM &,const SUP &);

      void constr_hess_I(double,const SUP &);

      static double gnorm(int);

      static int ghess2t(int,int);

      static int gt2hess(int,int);

      static void init();

      static void clear();

   private:

      //!list relating the single-particle space to the Hessian basis
      static vector< vector<int> > hess2t;

      //!list relating the single-particle space to the Hessian basis
      static int **t2hess;

      static vector<double> norm;

      //!hessian matrix
      Matrix *H;

      //!input gradient, output delta
      double *x;

};

#endif
