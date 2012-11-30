#ifndef TPSPNSM_H
#define TPSPNSM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "RecMat.h"

/**
 * @author Brecht Verstichel
 * @date 23-11-2012\n\n
 * This class TPSPnsM is a class written for the singly-traced Hessian matrix object,
 * being a rectangular matrix on TP for the rows, and SP for the columns, it inherits alle the function from its mother 
 * RecMat, a rectangular matrix, some special member functions and two lists 
 * The difference with TPSPM is that there is no symmetry between the two indices of SP and TP, hence the ns
 */
class TPSPnsM : public RecMat {

   public:
      
      //constructor
      TPSPnsM();

      //copy constructor
      TPSPnsM(const TPSPnsM &);

      //destructor
      virtual ~TPSPnsM();

      using RecMat::operator=;

      using RecMat::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d,int e,int z) const;

      void reorder(const DPM &);

   private:

};

#endif
