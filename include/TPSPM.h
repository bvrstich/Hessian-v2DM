#ifndef TPSPM_H
#define TPSPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

class TPTPnsM;

#include "RecMat.h"

/**
 * @author Brecht Verstichel
 * @date 23-11-2012\n\n
 * This class TPSPM is a class written for the singly-traced Hessian matrix object,
 * being a rectangular matrix on TP for the rows, and SP for the columns, it inherits alle the function from its mother 
 * RecMat, a rectangular matrix, some special member functions and two lists 
 */
class TPSPM : public RecMat {

   public:
      
      //constructor
      TPSPM();

      //copy constructor
      TPSPM(const TPSPM &);

      //destructor
      virtual ~TPSPM();

      using RecMat::operator=;

      using RecMat::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d,int e,int z) const;

      void dirprodtrace(double,const TPM &);

      void bar(double,const Hessian &);

      void bar(double,const TPTPnsM &);

      static int gspmm2s(int,int);

      static int gs2spmm(int,int);

      static int gspmmdim();
      
      static void init();

      static void clear();

   private:

      //!list relating the single-particle space to the column indices
      static vector< vector<int> > spmm2s;

      //!list relating the single-particle space to the column indices
      static int **s2spmm;

};

#endif
