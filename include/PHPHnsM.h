#ifndef PHPHnsM_H
#define PHPHnsM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

class PHSPnsM;

/**
 * @author Brecht Verstichel
 * @date 23-11-2012\n\n
 * This class PHPHnsM is a class written for matrices on ph-matrix space, it inherits alle the function from its mother Matrix
 */
class PHPHnsM : public Matrix {

   public:
      
      //constructor
      PHPHnsM();

      //copy constructor
      PHPHnsM(const PHPHnsM &);

      //destructor
      virtual ~PHPHnsM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d,int e,int z,int t,int h) const;

      double operator()(int I,int J,int K,int L) const;

      void square(const PHSPnsM &);

   private:

};

#endif
