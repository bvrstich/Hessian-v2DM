#ifndef PHPHM_H
#define PHPHM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

/**
 * @author Brecht Verstichel
 * @date 23-11-2012\n\n
 * This class PHPHM is a class written for matrices on ph-matrix space, it inherits alle the function from its mother Matrix
 */
class PHPHM : public Matrix {

   public:
      
      //constructor
      PHPHM();

      //copy constructor
      PHPHM(const PHPHM &);

      //destructor
      virtual ~PHPHM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d,int e,int z,int t,int h) const;

      double operator()(int I,int J,int K,int L) const;

      static int gphmm2ph(int,int);

      static int gph2phmm(int,int);

      static int gphmmdim();

      static void init();

      static void clear();

   private:

      //!list relating the PHPHM indices to the PH indices
      static vector< vector<int> > phmm2ph;

      //!list relating the PHPHM indices to the PH indices
      static int **ph2phmm;

};

#endif
