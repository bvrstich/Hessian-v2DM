#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

int **TPM::s2t;
vector< vector<int> > TPM::t2s;

double TPM::Sa;
double TPM::Sb;
double TPM::Sc;

/**
 * initialize the static lists
 */
void TPM::init(){

   int N = Tools::gN();
   int M = Tools::gM();

   //allocatie van sp2tp
   s2t = new int * [Tools::gM()];
   s2t[0] = new int [Tools::gM()*Tools::gM()];

   for(int i = 1;i < Tools::gM();++i)
      s2t[i] = s2t[i - 1] + Tools::gM();

   vector<int> v(2);

   //initialisatie van de twee arrays
   int t = 0;

   for(int a = 0;a < Tools::gM();++a)
      for(int b = a + 1;b < Tools::gM();++b){

         v[0] = a;
         v[1] = b;

         t2s.push_back(v);

         s2t[a][b] = t;
         s2t[b][a] = s2t[a][b];

         ++t;

      }

   Sa = 1.0;
   Sb = 0.0;
   Sc = 0.0;

#ifdef __Q_CON
   Sa += 1.0;
   Sb += (4.0*N*N + 2.0*N - 4.0*N*M + M*M - M)/(N*N*(N - 1.0)*(N - 1.0));
   Sc += (2.0*N - M)/((N - 1.0)*(N - 1.0));
#endif

#ifdef __G_CON
   Sa += 4.0;
   Sc += (2.0*N - M - 2.0)/((N - 1.0)*(N - 1.0));
#endif

#ifdef __T1_CON
   Sa += M - 4.0;
   Sb += (M*M*M - 6.0*M*M*N -3.0*M*M + 12.0*M*N*N + 12.0*M*N + 2.0*M - 18.0*N*N - 6.0*N*N*N)/( 3.0*N*N*(N - 1.0)*(N - 1.0) );
   Sc -= (M*M + 2.0*N*N - 4.0*M*N - M + 8.0*N - 4.0)/( 2.0*(N - 1.0)*(N - 1.0) );
#endif

#ifdef __T2_CON
   Sa += 5.0*M - 8.0;
   Sb += 2.0/(N - 1.0);
   Sc += (2.0*N*N + (M - 2.0)*(4.0*N - 3.0) - M*M)/(2.0*(N - 1.0)*(N - 1.0));
#endif

}

/**
 * deallocate the static lists
 */
void TPM::clear(){

   delete [] s2t[0];
   delete [] s2t;

}

/**
 * standard constructor: constructs Matrix object of dimension M*(M - 1)/2 and
 */
TPM::TPM() : Matrix(t2s.size()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix tpm_c
 * @param tpm_c object that will be copied into this.
 */
TPM::TPM(const TPM &tpm_c) : Matrix(tpm_c){ }

/**
 * destructor
 */
TPM::~TPM(){ }

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:\n\n
 * TPM(a,b,c,d) = -TPM(b,a,c,d) = -TPM(a,b,d,c) = TPM(b,a,c,d)
 * @param a first sp index that forms the tp row index i together with b
 * @param b second sp index that forms the tp row index i together with a
 * @param c first sp index that forms the tp column index j together with d
 * @param d second sp index that forms the tp column index j together with c
 * @return the number on place TPM(i,j) with the right phase.
 */
double TPM::operator()(int a,int b,int c,int d) const{

   if( (a == b) || (c == d) )
      return 0;
   else{

      int i = s2t[a][b];
      int j = s2t[c][d];

      int phase = 1;

      if(a > b)
         phase *= -1;
      if(c > d)
         phase *= -1;

      return phase*(*this)(i,j);

   }

}

ostream &operator<<(ostream &output,const TPM &tpm_p){

   for(int i = 0;i < tpm_p.gn();++i)
      for(int j = 0;j < tpm_p.gn();++j){

         output << i << "\t" << j << "\t|\t" << tpm_p.t2s[i][0] << "\t" << tpm_p.t2s[i][1]

            << "\t" << tpm_p.t2s[j][0] << "\t" << tpm_p.t2s[j][1] << "\t" << tpm_p(i,j) << endl;

      }

   return output;

}

/**
 * construct the hubbard hamiltonian with on site repulsion U
 * @param U onsite repulsion term
 * @param option == 0 use periodic boundary conditions, == 1 use no pbc
 */
void TPM::hubbard(int option,double U){

   int a,b,c,d;//sp orbitals

   double ward = 1.0/(Tools::gN() - 1.0);

   for(int i = 0;i < gn();++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < gn();++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = 0;

         if(option == 0){//pbc

            //eerst hopping
            if( (a == c) && ( ( (b + 2)%Tools::gM() == d ) || ( b == (d + 2)%Tools::gM() ) ) )
               (*this)(i,j) -= ward;

            if( (b == c) && ( ( (a + 2)%Tools::gM() == d ) || ( a == (d + 2)%Tools::gM() ) ) )
               (*this)(i,j) += ward;

            if( (b == d) && ( ( (a + 2)%Tools::gM() == c ) || ( a == (c + 2)%Tools::gM() ) ) )
               (*this)(i,j) -= ward;

         }
         else{//no pbc

            //eerst hopping
            if( (a == c) && ( ( (b + 2) == d ) || ( b == (d + 2) ) ) )
               (*this)(i,j) -= ward;

            if( (b == c) && ( ( (a + 2) == d ) || ( a == (d + 2) ) ) )
               (*this)(i,j) += ward;

            if( (b == d) && ( ( (a + 2) == c ) || ( a == (c + 2) ) ) )
               (*this)(i,j) -= ward;

         }

         //on site interaction
         if( (a % 2) == 0 && (c % 2) == 0 )
            if(a == (b - 1) && c == (d - 1) && a == c)
               (*this)(i,j) += U;

      }

   }

   this->symmetrize();

}

/**
 * The Q map
 * @param option = 1, regular Q map , = -1 inverse Q map
 * @param tpm_d the TPM of which the Q map is taken and saved in this.
 */
void TPM::Q(int option,const TPM &tpm_d){

   double a = 1;
   double b = 1.0/(Tools::gN()*(Tools::gN() - 1.0));
   double c = 1.0/(Tools::gN() - 1.0);

   this->Q(option,a,b,c,tpm_d);

}

/**
 * The Q-like map: see primal-dual.pdf for more info (form: Q(A,B,C)(TPM) )
 * @param option = 1, regular Q-like map , = -1 inverse Q-like map
 * @param A factor in front of the two particle piece of the map
 * @param B factor in front of the no particle piece of the map
 * @param C factor in front of the single particle piece of the map
 * @param tpm_d the TPM of which the Q-like map is taken and saved in this.
 */
void TPM::Q(int option,double A,double B,double C,const TPM &tpm_d){

   if(option == -1){

      B = (B*A + B*C*Tools::gM() - 2.0*C*C)/( A * (C*(Tools::gM() - 2.0) -  A) * ( A + B*Tools::gM()*(Tools::gM() - 1.0) - 2.0*C*(Tools::gM() - 1.0) ) );
      C = C/(A*(C*(Tools::gM() - 2.0) - A));
      A = 1.0/A;

   }

   SPM spm;

   //de trace*2 omdat mijn definitie van trace in berekeningen over alle (alpha,beta) loopt
   double ward = B*tpm_d.trace()*2.0;

   //construct de spm met schaling C
   spm.bar(C,tpm_d);

   for(int i = 0;i < gn();++i){

      int a = t2s[i][0];
      int b = t2s[i][1];

      for(int j = i;j < gn();++j){

         int c = t2s[j][0];
         int d = t2s[j][1];

         (*this)(i,j) = A*tpm_d(i,j);

         if(i == j)
            (*this)(i,i) += ward;

         if(a == c)
            (*this)(i,j) -= spm(b,d);

         if(b == c)
            (*this)(i,j) += spm(a,d);

         if(b == d)
            (*this)(i,j) -= spm(a,c);

      }
   }

   this->symmetrize();

}

/**
 * initialize this onto the unitmatrix with trace N*(N - 1)/2
 */
void TPM::unit(){

   double ward = Tools::gN()*(Tools::gN() - 1.0)/(Tools::gM() * (Tools::gM() - 1.0));

   for(int i = 0;i < gn();++i){

      (*this)(i,i) = ward;

      for(int j = i + 1;j < gn();++j)
         (*this)(i,j) = (*this)(j,i) = 0.0;

   }

}

/**
 * orthogonal projection onto the space of traceless matrices
 */
void TPM::proj_Tr(){

   double ward = 2.0 * (this->trace())/(double)(Tools::gM() * (Tools::gM() - 1.0));

   for(int i = 0;i < gn();++i)
      (*this)(i,i) -= ward;

}

/**
 * Deduct the unitmatrix times a constant (scale) from this.\n\n
 * this -= scale* 1
 * @param scale the constant
 */
void TPM::min_unit(double scale){

   for(int i = 0;i < gn();++i)
      (*this)(i,i) -= scale;

}

void TPM::in_sp(const char *filename){

   ifstream input(filename);

   double value;

   int a,b,c,d;

   int i,j;

   while(input >> a >> b >> c >> d >> value){

      i = s2t[a][b];
      j = s2t[c][d];

      (*this)(i,j) = value;

   }

   this->symmetrize();

}

/**
 * @return the expectation value of the size of the spin: S^2
 */
double TPM::S_2() const{

   //first diagonal elements:
   int a,b;
   double s_a,s_b;

   double ward = 0.0;

   for(int i = 0;i < gn();++i){

      a = t2s[i][0];
      b = t2s[i][1];

      s_a = ( 1.0 - 2 * (a % 2) )/2;
      s_b = ( 1.0 - 2 * (b % 2) )/2;

      ward += ( (1 + s_a*s_a + s_b*s_b)/(Tools::gN() - 1.0) + 2*s_a*s_b ) * (*this)(i,i);

   }

   //then the off diagonal elements: a and b are sp indices
   for(int a = 0;a < Tools::gM()/2;++a)
      for(int b = 0;b < Tools::gM()/2;++b)
         ward += (*this)(2*a,2*b + 1,2*a + 1,2*b);

   return ward;

}

/**
 * fill the TPM object with the S^2 matrix
 */
void TPM::set_S_2(){

   int a,b,c,d;

   double s_a,s_b;

   for(int i = 0;i < gn();++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < gn();++j){

         c = t2s[j][0];
         d = t2s[j][1];

         //init
         (*this)(i,j) = 0.0;

         if(i == j){//diagonal stuff

            s_a = ( 1.0 - 2 * (a % 2) )/2;
            s_b = ( 1.0 - 2 * (b % 2) )/2;

            (*this)(i,i) = (1 + s_a*s_a + s_b*s_b)/(Tools::gN() - 1.0) + 2*s_a*s_b;

            if(a/2 == b/2 && a % 2 == 0 && b % 2 == 1)
               (*this)(i,i) -= 1.0;

         }

         //then the off-diagonal elements
         if(a % 2 == 0 && b % 2 == 1 && a/2 != b/2)//a up and b down
            if(a + 1 == c && b == d + 1)
               (*this)(i,j) += 1.0;

      }

   }

   this->symmetrize();

}

/**
 * @return the dimension of a TPM matrix
 */
int TPM::gn(){

   return t2s.size();

}

/**
 * convert the solution Newton object to a TPM Matrix
 * @param newton input Newton object
 */
void TPM::convert(const Newton &newton){

   for(int i = 0;i < gn();++i)
      for(int j = i;j < gn();++j)
         (*this)(i,j) = newton.gx(i,j)/(2.0*Newton::gnorm(i,j));

   this->symmetrize();

}

/**
 * access to the lists from outside the class
 */
int TPM::gt2s(int i,int option){

   return t2s[i][option];

}

/**
 * access to the lists from outside the class
 */
int TPM::gs2t(int a,int b){

   return s2t[a][b];

}

/**
 * The G-map that maps a PHM object onto a TPM object.
 * @param option = 1, G_down - map is used, = -1 G^{-1}_up - map is used.
 * @param phm input PHM 
 */
void TPM::G(int option,const PHM &phm){

   SPM spm;

   if(option == 1)
      spm.bar(1.0/(Tools::gN() - 1.0),phm);
   else
      spm.bar(1.0/(Tools::gM() - Tools::gN() + 1.0),phm);

   int a,b,c,d;

   for(int i = 0;i < gn();++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < gn();++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = phm(b,d,c,a) - phm(a,d,c,b) - phm(b,c,d,a) + phm(a,c,d,b);

         if(b == d)
            (*this)(i,j) += spm(a,c);

         if(b == c)
            (*this)(i,j) -= spm(a,d);

         if(a == c)
            (*this)(i,j) += spm(b,d);

      }

   }

   //nog schalen met 4 voor G^{-1}
   if(option == -1)
      this->dscal(0.25);

   this->symmetrize();

}

/**
 * map a DPM (dpm) on a TPM (*this) with a T1 map, (Q-like map), watch out for the inverse
 * up map, when M = 2*N it is singular! So don't use it!:
 * @param option = +1 T1_down , =-1 inverse T1_up
 * @param dpm The input DPM
 */
void TPM::T(int option,const DPM &dpm){

   TPM tpm;
   tpm.bar(1.0,dpm);

   if(option == 1){

      double a = 1;
      double b = 1.0/(3.0*Tools::gN()*(Tools::gN() - 1.0));
      double c = 0.5/(Tools::gN() - 1.0);

      this->Q(1,a,b,c,tpm);

   }
   else{//option == -1

      double a = Tools::gM() - 4.0;
      double b = (Tools::gM() - Tools::gN() - 2.0)/(Tools::gN()*(Tools::gN() - 1.0));
      double c = (Tools::gM() - Tools::gN() - 2.0)/(Tools::gN() - 1.0);

      this->Q(-1,a,b,c,tpm);

   }

}

/**
 * calculate the trace of one pair of sp indices of a DPM an put in (*this):\n\n
 * TPM(a,b,d,e) = sum_{c} DPM(a,b,c,d,e,c)
 * @param dpm input DPM
 */
void TPM::bar(double scale,const DPM &dpm){

   int a,b,c,d;

   for(int i = 0;i < gn();++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < gn();++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = 0.0;

         for(int l = 0;l < Tools::gM();++l)
            (*this)(i,j) += dpm(a,b,l,c,d,l);

         (*this)(i,j) *= scale;

      }
   }

   this->symmetrize();

}

/**
 * Map a PPHM (pphm) onto a TPM object (*this) with a T2 down map, see primal_dual.pdf for more information
 * @param pphm input PPHM
 */
void TPM::T(const PPHM &pphm){

   int M = Tools::gM();
   int N = Tools::gN();

   //first make some necessary derivate matrices of pphm
   TPM bar;
   bar.bar(pphm);

   PHM phm;
   phm.bar(pphm);

   //watch out, scaling for spm is not the usual!
   SPM spm;

   for(int a = 0;a < M;++a)
      for(int b = a;b < M;++b){

         spm(a,b) = 0;

         for(int c = 0;c < M;++c)
            spm(a,b) += phm(c,a,c,b);

         spm(a,b) *= 0.5/(N - 1.0);

      }

   int a,b,c,d;

   for(int i = 0;i < gn();++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < gn();++j){

         c = t2s[j][0];
         d = t2s[j][1];

         //first the tp part:
         (*this)(i,j) = bar(i,j);

         //then the ph part:
         (*this)(i,j) -= phm(d,a,b,c) - phm(d,b,a,c) - phm(c,a,b,d) + phm(c,b,a,d);

         //finaly the three sp parts:
         if(b == d)
            (*this)(i,j) += spm(a,c);

         if(b == c)
            (*this)(i,j) -= spm(a,d);

         if(a == c)
            (*this)(i,j) += spm(b,d);

      }
   }

   this->symmetrize();

}

/**
 * Map a PPHM (pphm) object on a TPM (*this) object by tracing one pair of indices from the pphm (for more info, see primal_dual.pdf)
 * @param pphm input PPHM
 */
void TPM::bar(const PPHM &pphm){

   int M = Tools::gM();

   int a,b,c,d;

   for(int i = 0;i < gn();++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < gn();++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = 0.0;

         for(int l = 0;l < M;++l)
            (*this)(i,j) += pphm(a,b,l,c,d,l);

      }
   }

   this->symmetrize();

}

/**
 * Collaps a SUP matrix S onto a TPM matrix like this:\n\n
 * sum_i Tr (S u^i)f^i = this
 * @param option = 0, project onto full symmetric matrix space, = 1 project onto traceless symmetric matrix space
 * @param S input SUP
 */
void TPM::collaps(int option,const SUP &S){

   *this = S.gI();

   TPM hulp;

   hulp.Q(1,S.gQ());

   *this += hulp;

#ifdef __G_CON
   hulp.G(1,S.gG());

   *this += hulp;
#endif

#ifdef __T1_CON
   hulp.T(1,S.gT1());

   *this += hulp;
#endif

#ifdef __T2_CON
   hulp.T(S.gT2());

   *this += hulp;
#endif

   if(option == 1)
      this->proj_Tr();

}

/**
 * ( Overlapmatrix of the U-basis ) - map, maps a TPM onto a different TPM, this map is actually a Q-like map
 * for which the paramaters a,b and c are calculated in primal_dual.pdf. Since it is a Q-like map the inverse
 * can be taken as well.
 * @param option = 1 direct overlapmatrix-map is used , = -1 inverse overlapmatrix map is used
 * @param tpm_d the input TPM
 */
void TPM::S(int option,const TPM &tpm_d){

   this->Q(option,Sa,Sb,Sc,tpm_d);

}
