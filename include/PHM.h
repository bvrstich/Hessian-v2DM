#ifndef PHM_H
#define PHM_H

#include <iostream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"
#include "TPM.h"

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class, PHM, is a class written for particle-hole matrices, it inherits all the functions from its mother class
 * Matrix, some special member functions and two lists that give the relationship between the sp and the ph basis.
 */
class PHM : public Matrix {

    /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << phm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << phm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param phm_p the TPM you want to print
    */
   friend ostream &operator<<(ostream &output,PHM &phm_p);

   public:
      
      //constructor
      PHM();

      //copy constructor
      PHM(const PHM &);

      //destructor
      virtual ~PHM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d) const;

      //change the numbers in sp mode
      double &operator()(int a,int b,int c,int d);

      void G(int option,const TPM &);

      static int gs2ph(int,int);

      static int gph2s(int,int);

      static int gn();

      static void init();

      static void clear();

   private:

      //!static list of dimension [n_ph][2] that takes in a ph index i and returns two sp indices: a = ph2s[i][0] and b = ph2s[i][1]
      static vector< vector<int> > ph2s;

      //!static list of dimension [M][M] that takes two sp indices a,b and returns a ph index i: i = s2t[a][b]
      static int **s2ph;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/