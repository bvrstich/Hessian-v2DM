#ifndef PHSPNSM_H
#define PHSPNSM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "RecMat.h"

/**
 * @author Brecht Verstichel
 * @date 23-11-2012\n\n
 * This class PHSPnsM is a rectangular matrix on PH for the rows, and SP for the columns, it inherits alle the function from its mother 
 * RecMat, a rectangular matrix, some special member functions and two lists 
 * The difference with PHSPM is that there is no symmetry between the two indices of SP and PH, hence the ns
 */
class PHSPnsM : public RecMat {

   public:
      
      //constructor
      PHSPnsM();

      //copy constructor
      PHSPnsM(const PHSPnsM &);

      //destructor
      virtual ~PHSPnsM();

      using RecMat::operator=;

      using RecMat::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d,int e,int z) const;

      void reorder(const PPHM &);

      static int gphmm2ph(int,int);

      static int gph2phmm(int,int);

      static int gphmmdim();

      static void init();

      static void clear();

   private:

      //!list relating the non-symmetrical PHPHM indices to the PH indices
      static vector< vector<int> > phmm2ph;

      //!list relating the non-symmetrical PHPHM indices to the PH indices
      static int **ph2phmm;

};

#endif
