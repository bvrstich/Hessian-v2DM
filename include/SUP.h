#ifndef SUP_H
#define SUP_H

#include <iostream> 
#include <fstream> 

using std::ostream;
using std::ofstream;
using std::ifstream;

#include "TPM.h"

#ifdef PQ

#define __Q_CON

#endif

#ifdef PQG

#define __Q_CON
#define __G_CON

#endif

#ifdef PQGT1

#define __Q_CON
#define __G_CON
#define __T1_CON

#endif

#ifdef PQGT2

#define __Q_CON
#define __G_CON
#define __T2_CON

#endif

#ifdef PQGT

#define __Q_CON
#define __G_CON
#define __T1_CON
#define __T2_CON

#endif

class EIG;

/**
 * @author Brecht Verstichel
 * @date 09-03-2010\n\n
 * This class, SUP is a blockmatrix over the carrierspace's of active N-representability conditions. 
 * This class contains two TPM objects, and if compiled with the right option a PHM, DPM or PPHM object, 
 */
class SUP{
  
   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << sup_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << sup_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param sup_p the SUP you want to print
    */
   friend ostream &operator<<(ostream &output,const SUP &sup_p);

   public:

      //constructor
      SUP();

      //copy constructor
      SUP(const SUP &);

      //destructor
      ~SUP();

      //overload += operator
      SUP &operator+=(const SUP &);

      //overload -= operator
      SUP &operator-=(const SUP &);

      //overload equality operator
      SUP &operator=(const SUP &);

      //overload equality operator
      SUP &operator=(double &);

      double ddot(const SUP &) const;

      void invert();

      void dscal(double alpha);

      void sqrt(int option);

      void L_map(const SUP &,const SUP &);

      void daxpy(double alpha,const SUP &);

      SUP &mprod(const SUP &,const SUP &);

      void fill(const TPM &);

      void fill();

      void fill_Random();

      TPM &gI();

      const TPM &gI() const;

#ifdef __Q_CON
      TPM &gQ();

      const TPM &gQ() const;
#endif

#ifdef __G_CON
      PHM &gG();

      const PHM &gG() const;
#endif

#ifdef __T1_CON
      DPM &gT1();

      const DPM &gT1() const;
#endif

#ifdef __T2_CON
      PPHM &gT2();

      const PPHM &gT2() const;
#endif

   private:

      //!pointer to the TPM object containing the I contribution to the SUP matrix
      TPM *I;

#ifdef __Q_CON
      //!pointer to the TPM object containing the Q contribution to the SUP matrix
      TPM *Q;
#endif

#ifdef __G_CON
      //!pointer to the particle hole matrix
      PHM *G;
#endif

#ifdef __T1_CON
      //!pointer to the three-particle matrix DPM
      DPM *T1;
#endif

#ifdef __T2_CON
      //!pointer to the two-particle-one-hole matrix PPHM
      PPHM *T2;
#endif

};

#endif
