#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "include.h"

/**
 * constructor, makes matrix of dimension M
 */
SPM::SPM() : Matrix(Tools::gM()) { }

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
SPM::SPM(const SPM &spm_copy) : Matrix(spm_copy) { }

/**
 * destructor
 */
SPM::~SPM(){ }

ostream &operator<<(ostream &output,const SPM &spm_p){

   for(int i = 0;i < Tools::gM();++i)
      for(int j = 0;j < Tools::gM();++j)
         output << i << "\t" << j << "\t" << spm_p(i,j) << endl;

   return output;

}

/**
 * construct the SPM object by tracing out one pair of indices from a TPM object
 * @param tpm input TPM object
 */
void SPM::bar(double scale,const TPM &tpm){

   for(int a = 0;a < Tools::gM();++a)
      for(int b = a;b < Tools::gM();++b){

         (*this)(a,b) = 0.0;

         for(int c = 0;c < Tools::gM();++c)
            (*this)(a,b) += tpm(a,c,b,c);

         (*this)(a,b) *= scale;

      }
   
   this->symmetrize();

}
