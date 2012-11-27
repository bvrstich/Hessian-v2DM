#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ifstream;
using std::endl;

#include "include.h"

vector< vector<int> > PHM::ph2s;
int **PHM::s2ph;

/**
 * initialize the static variables and allocate the static lists
 */
void PHM::init(){

   //allocate
   s2ph = new int * [Tools::gM()];

   for(int a = 0;a < Tools::gM();++a)
      s2ph[a] = new int [Tools::gM()];

   int ph = 0;

   vector<int> v(2);

   for(int a = 0;a < Tools::gM();++a)
      for(int b = 0;b < Tools::gM();++b){

         s2ph[a][b] = ph;

         v[0] = a;
         v[1] = b;

         ph2s.push_back(v);

         ++ph;

      }

}

/**
 * deallocate the lists
 */
void PHM::clear(){

   for(int a = 0;a < Tools::gM();++a)
      delete [] s2ph[a];

   delete [] s2ph;

}

/**
 * standard constructor: constructs Matrix object of dimension M*M and
 */
PHM::PHM() : Matrix(ph2s.size()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*M and copies the content of phm_c into it,
 * @param phm_c PHM to be copied into (*this)
 */
PHM::PHM(const PHM &phm_c) : Matrix(phm_c){ }

/**
 * destructor
 */
PHM::~PHM(){ }

/**
 * access the elements of the matrix in sp mode, 
 * @param a first sp index that forms the ph row index i together with b
 * @param b second sp index that forms the ph row index i together with a
 * @param c first sp index that forms the ph column index j together with d
 * @param d second sp index that forms the ph column index j together with c
 * @return the number on place PHM(i,j)
 */
double &PHM::operator()(int a,int b,int c,int d){

   int i = s2ph[a][b];
   int j = s2ph[c][d];

   return (*this)(i,j);

}

/**
 * access the elements of the matrix in sp mode, 
 * @param a first sp index that forms the ph row index i together with b
 * @param b second sp index that forms the ph row index i together with a
 * @param c first sp index that forms the ph column index j together with d
 * @param d second sp index that forms the ph column index j together with c
 * @return the number on place PHM(i,j)
 */
double PHM::operator()(int a,int b,int c,int d) const
{
   int i = s2ph[a][b];
   int j = s2ph[c][d];

   return (*this)(i,j);
}

ostream &operator<<(ostream &output,PHM &phm_p){

   for(int i = 0;i < phm_p.gn();++i)
      for(int j = 0;j < phm_p.gn();++j){

         output << i << "\t" << j << "\t|\t" << phm_p.ph2s[i][0] << "\t" << phm_p.ph2s[i][1]

            << "\t" << phm_p.ph2s[j][0] << "\t" << phm_p.ph2s[j][1] << "\t" << phm_p(i,j) << endl;

      }

   return output;

}

/**
 * De G map, maps a TPM object on a PHM object.
 * @param option = 1 G_up map is used, = -1 G^{-1}_down map is used
 * @param tpm input TPM
 */
void PHM::G(int option,const TPM &tpm)
{
   SPM spm;

   if(option == 1)
      spm.bar(1.0/(Tools::gN() - 1.0),tpm);
   else
      spm.bar(1.0/(Tools::gM() - Tools::gN() + 1.0),tpm);

   int a,b,c,d;

   for(int i = 0;i < gn();++i){

      a = ph2s[i][0];
      b = ph2s[i][1];

      for(int j = i;j < gn();++j){

         c = ph2s[j][0];
         d = ph2s[j][1];

         (*this)(i,j) = -tpm(a,d,c,b);

         if(b == d)
            (*this)(i,j) += spm(a,c);

      }
   }
   
   //nog schalen met 4 voor G^{-1}, door de G down die eigenlijk een factor 4 te groot is
   if(option == -1)
      this->dscal(0.25);

   this->symmetrize();

}

/* vim: set ts=3 sw=3 expandtab :*/
