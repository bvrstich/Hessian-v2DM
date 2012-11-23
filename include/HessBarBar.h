#ifndef HESSBARBAR_H
#define HESSBARBAR_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

/**
 * @author Brecht Verstichel
 * @date 23-11-2012\n\n
 * This class HessBarBar is a class written for the singly-traced Hessian matrix object, it inherits alle the function from its mother 
 * Matrix, some special member functions and two lists 
 */
class HessBarBar : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param hessbb_p the HessBarBar you want to print
    */
   friend ostream &operator<<(ostream &output,const HessBarBar &hessbb_p);

   public:
      
      //constructor
      HessBarBar();

      //copy constructor
      HessBarBar(const HessBarBar &);

      //destructor
      virtual ~HessBarBar();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d) const;

      void dirprodtrace(double,const TPM &);

      static int ghb2s(int,int);

      static int gs2hb(int,int,int);
      
      static void init();

      static void clear();

   private:

      //!list relating the single-particle space to the HessBarBar basis
      static vector< vector<int> > hb2s;

      //!list relating the single-particle space to the HessBarBar basis
      static int ***s2hb;

};

#endif
