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

      //destructor
      virtual ~Newton();

      const Hessian &gH() const;

      const double *gx() const;

      double *gx();

      double gx(int,int) const;

      double gx(int,int,int,int) const;

      void construct(double,const TPM &,const SUP &);

      void gradient(double,const TPM &,const SUP &);

      static double gnorm(int);

      static double gnorm(int,int);

      static void init();

      static void clear();

   private:

      static double *norm;

      //!hessian matrix
      Hessian *H;

      //!input gradient, output delta
      double *x;

};

#endif
