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

/**
 * initialize the static lists
 */
void TPM::init(){

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
 * perform a line search what step size in along the Newton direction is ideal.
 * @param t potential scaling factor
 * @param P SUP matrix containing the inverse of the constraints (carrier space matrices)
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double TPM::line_search(double t,SUP &P,const TPM &ham){

   double tolerance = 1.0e-5*t;

   if(tolerance < 1.0e-12)
      tolerance = 1.0e-12;

   //neem de wortel uit P
   P.sqrt(1);

   //maak eerst een SUP van delta
   SUP S_delta;

   S_delta.fill(*this);

   //hulpje om dingskes in te steken:
   SUP hulp;

   hulp.L_map(P,S_delta);

   EIG eigen(hulp);

   double a = 0;

   double b = -1.0/eigen.min();

   double c(0);

   double ham_delta = ham.ddot(*this);

   while(b - a > tolerance){

      c = (b + a)/2.0;

      if( (ham_delta - t*eigen.lsfunc(c)) < 0.0)
         a = c;
      else
         b = c;

   }

   return c;

}

/**
 * perform a line search what step size in along the Newton direction is ideal, this one is used for extrapolation.
 * @param t potential scaling factor
 * @param rdm TPM containing the current approximation of the rdm.
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double TPM::line_search(double t,const TPM &rdm,const TPM &ham){

   SUP P;

   P.fill(rdm);

   P.invert();

   return this->line_search(t,P,ham);

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
