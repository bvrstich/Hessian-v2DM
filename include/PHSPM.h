#ifndef PHSPM_H
#define PHSPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "RecMat.h"

class PHM;

/**
 * @author Brecht Verstichel
 * @date 23-11-2012\n\n
 * This class PHSPM is a class written for the singly-traced symmetrized direct product of ph-matrices
 * being a rectangular matrix on PH for the rows, and SP for the columns, it inherits alle the function from its mother 
 * RecMat, a rectangular matrix, some special member functions and two lists 
 */
class PHSPM : public RecMat {

   public:
      
      //constructor
      PHSPM();

      //copy constructor
      PHSPM(const PHSPM &);

      //destructor
      virtual ~PHSPM();

      using RecMat::operator=;

      using RecMat::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d,int e,int z) const;

      void dirprodtrace(double,const PHM &);

      static int gphmm2ph(int,int);

      static int gph2phmm(int,int);

      static int gphmmdim();
      
      static void init();

      static void clear();

   private:

      //!list relating the ph-matrix space indices to the ph-indices (like hess2t)
      static vector< vector<int> > phmm2ph;

      //!list relating the single-particle space to the column indices
      static int **ph2phmm;

};

#endif
