#ifndef TPTPnsM_H
#define TPTPnsM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

class TPSPnsM;

/**
 * @author Brecht Verstichel
 * @date 23-11-2012\n\n
 * This class TPTPnsM is a class written for matrices on tp-matrix space, it inherits alle the function from its mother Matrix
 */
class TPTPnsM : public Matrix {

   public:
      
      //constructor
      TPTPnsM();

      //copy constructor
      TPTPnsM(const TPTPnsM &);

      //destructor
      virtual ~TPTPnsM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d,int e,int z,int t,int h) const;

      double operator()(int I,int J,int K,int L) const;

      void square(const TPSPnsM &);

      static int gtpmm2t(int,int);

      static int gt2tpmm(int,int);

      static void init();

      static void clear();

      static int gn();

   private:

      //!list relating the single-particle space to the Hessian basis
      static vector< vector<int> > tpmm2t;

      //!list relating the single-particle space to the Hessian basis
      static int **t2tpmm;

};

#endif
