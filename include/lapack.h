#ifndef LAPACK_H
#define LAPACK_H

//lapack definities
extern "C" {

   void dgeqrf_(int *m,int *n,double *A,int *LDA,double *tau,double *WORK,int *LWORK,int *INFO);
   void dorgqr_(int *m,int *n, int *k,double *A,int  *LDA,double *tau,double *WORK,int *LWORK,int *INFO);
   void dgelqf_(int *m,int *n,double *A,int *LDA,double *tau,double *WORK,int *LWORK,int *INFO);
   void dorglq_(int *m,int *n, int *k,double *A,int  *LDA,double *tau,double *WORK,int *LWORK,int *INFO);
   void dcopy_(int *n,double *x,int *incx,double *y,int *incy);
   void daxpy_(int *n,double *alpha,double *x,int *incx,double *y,int *incy);
   void dscal_(int *n,const double *alpha,double *x,int *incx);
   void dgemm_(char *transA,char *transB,int *m,int *n,int *k,double *alpha,double *A,int *lda,double *B,int *ldb,double *beta,double *C,int *ldc);
   void dsymm_(char *side,char *uplo,int *m,int *n,double *alpha,double *A,int *lda,double *B,int *ldb,double *beta,double *C,int *ldc);
   void dgemv_(char *trans,int *m,int *n,double *alpha,double *A,int *lda,double *x,int *incx,double *beta,double *y,int *incy);
   double ddot_(int *n,double *x,int *incx,double *y,int *incy);
   void dsyev_(char *jobz,char *uplo,int *n,double *A,int *lda,double *W,double *work,int *lwork,int *info);
   void dpotrf_(char *uplo,int *n,double *A,int *lda,int *INFO);
   void dpotri_(char *uplo,int *n,double *A,int *lda,int *INFO);
   void dgetrf_(int *m ,int *n,double *A,int *lda,int *ipiv,int *INFO);
   void dgetrs_(char *trans ,int *n,int *nrhs,double *A,int *lda,int *ipiv,double *B,int *ldb,int *INFO);
   void dsytrf_(char *uplo,int *n,double *A,int *lda,int *ipiv,double *work,int *lwork,int *INFO);
   void dsytrs_(char *uplo,int *n,int *nrhs,double *A,int *lda,int *ipiv,double *B,int *ldb,int *INFO);
   void dpotrs_(char *uplo,int *n,int *nrhs,double *A,int *lda,double *B,int *ldb,int *INFO);
   void dsymv_(char *uplo,int *n,double *alpha,double *A,int *lda,const double *x,int *incx,double *beta,double *y,int *incy);

}

#endif
