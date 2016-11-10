#ifndef _MATRIX_OPERATION_H_
#define _MATRIX_OPERATION_H_
//#include "nb.h"
/*******************   eigen problem and SVD   *****************/
void exponential_symmatrix(DoubleMatrix &exp_M, const DoubleMatrix &M, double tau);
void exponential_symmatrix(ComplexMatrix &exp_M, const ComplexMatrix &M, std::complex<double> t);
void eigen_symmetric_solver(const DoubleMatrix &h, DoubleVector &eigval, DoubleMatrix &eigvec, char jobz='N', int ill=1); 
void eigen_symmetric_solver(const ComplexMatrix &h, DoubleVector &eigval, ComplexMatrix &eigvec, char jobz='N', int ill=1);
void eigen_general_solver(const DoubleMatrix &h, DoubleVector &lambda, DoubleMatrix &U_in, DoubleMatrix &V_in, char lOrR);
void eigen_general_solver(const ComplexMatrix &h, ComplexVector &lambda, ComplexMatrix &U_in, ComplexMatrix &V_in, char LOrR);
void svd_solver(const ComplexMatrix &h_in, ComplexMatrix &u_out, DoubleVector &s_out,ComplexMatrix &vt_out);
void svd_solver(const DoubleMatrix &h_in, DoubleMatrix &u_out, DoubleVector &s_out,DoubleMatrix &vt_out);

void exponential_symmatrix(DoubleMatrix &exp_M, const DoubleMatrix &M, double tau){
  /*
   * exp_M = exp(-tau M);
   */
  int order=M.size();
  DoubleVector eigval(order);
  DoubleMatrix eigvec(order);
  for(int i=0; i<order; i++)
    eigvec[i].resize(order);
  eigen_symmetric_solver(M,eigval,eigvec,'V');
  for(int i=0; i < order; i++)
    for(int j=0; j < order; j++){
      exp_M[i][j]=0;
      for(int k=0; k < order; k++)
        exp_M[i][j] += exp(-eigval[k]*tau)*eigvec[i][k]*eigvec[j][k];
    }
}
void exponential_symmatrix(DoubleMatrix &M, double tau){
  /*
   * M = exp(-tau M)
   */
  exponential_symmatrix(M,M,tau);
}

void exponential_symmatrix(ComplexMatrix &exp_M, const ComplexMatrix &M, std::complex<double> t){
  /*
   * exp_M = exp(-i*t*A);
   */
  int order=M.size();
  DoubleVector eigval(order);
  ComplexMatrix eigvec(order);
  for(int i=0; i < order; i++)
    eigvec[i].resize(order);
  eigen_symmetric_solver(M,eigval,eigvec,'V');
  std::complex<double> eye(0,1.0);
  for(int i=0; i < order; i++)
    for(int j=0; j < order; j++){
      exp_M[i][j]=0;
      for(int k=0; k < order; k++)
        exp_M[i][j] += exp(-eye*eigval[k]*t)*eigvec[i][k]*std::conj(eigvec[j][k]);
    }
}
void exponential_symmatrix(ComplexMatrix &M, std::complex<double> t){
  exponential_symmatrix(M,M,t);
}

void eigen_symmetric_solver(const DoubleMatrix &h, DoubleVector &eigval, DoubleMatrix &eigvec, char jobz, int ill){
  /**
   *[in]
   *jobz: 'N':  Compute eigenvalues only;
   *      'V':  Compute eigenvalues and eigenvectors.
   *[out]
   * m: total number of eigenvalues found, iu-il+1
   * w: DOUBLE PRECISION array, dimension (N)
   *    On normal exit, the first M elements contain the selected eigenvalues in ascending order.
   * z: the first M columns give the selected eigenvalues of matrix A
  **/
  if(ill<1){
    std::cerr << "eigen_symmetric_solver: il >= 1";
    std::exit(1.0);
  }
  lapack_int n, lda, il, iu, m, ldz, info;
  double vl, vu, abstol;
  double *a, *w, *z;
  int *ifail;  
  n = h.size();lda=n; il=ill; iu=eigval.size()+il-1; abstol=0; ldz=n;
  
  a = (double *) calloc(n*n,sizeof(double));
  w = (double *) calloc(n,sizeof(double));
  z = (double *) calloc(ldz*(iu-il+1),sizeof(double));
  ifail = (int *) calloc(n,sizeof(int));
  
  for(int j=0; j<n; j++)
    for(int i=0; i<n; i++)
         a[i+j*n] = h[i][j]; // transform matrix h to col_major

  info = LAPACKE_dsyevx(LAPACK_COL_MAJOR,jobz,'I','U',n,a,lda,vl,vu,il,iu,abstol,&m,w,z,ldz,ifail);
  if(info == 0){
    for(int i=0; i<iu-il+1; i++)
      eigval[i]= w[i];
    if(jobz=='V')
      for(int i=0; i<iu-il+1; i++)
        for(int j=0; j<n; j++)
          eigvec[j][i]= z[j+i*n];
  }else{
    std::cerr << -info-1 << "[double] LAPACKE_dsyevs: eigval is not converge" << std::endl;
  }
  free(a);
  free(w);
  free(z);
  free(ifail);
}

void eigen_symmetric_solver(const ComplexMatrix &h, DoubleVector &eigval, ComplexMatrix &eigvec, char jobz, int ill){
  /**
   * ZHEEVX computes selected eigenvalues and, optionally, eigenvectors
   * of a complex Hermitian matrix A.  Eigenvalues and eigenvectors can
   * be selected by specifying either a range of values or a range of
   * indices for the desired eigenvalues.
  **/
  if(ill<1){
    std::cerr << "eigen_symmetric_solver: il >= 1";
    std::exit(1.0);
  }
  lapack_int n,lda,il,iu,m,ldz,info;
  double vl,vu,abstol;
  double *w;
  lapack_complex_double *a,*z;
  lapack_int *ifail;
  n = h.size(); il=ill; iu=eigval.size()+il-1; lda=n; abstol=0; ldz=n;  
  
  a = (lapack_complex_double *) calloc(n*n,sizeof(lapack_complex_double));
  w = (double *) calloc(n,sizeof(double));
  z = (lapack_complex_double *) calloc(ldz*(iu-il+1),sizeof(lapack_complex_double));
  ifail = (int *) calloc(n,sizeof(int));
  
  for(int j=0; j<n; j++)
    for(int i=0; i<n; i++)
      a[i+j*n] = lapack_make_complex_double(h[i][j].real(),h[i][j].imag()); // transform matrix h to col_major

  info = LAPACKE_zheevx(LAPACK_COL_MAJOR,jobz,'I','U',n,a,lda,vl,vu,il,iu,abstol,&m,w,z,ldz,ifail);
  if(info == 0){
    for(int i=0; i<iu-il+1; i++)
      eigval[i]= w[i];
    if(jobz=='V')
      for(int i=0; i<iu-il+1; i++)
        for(int j=0; j<n; j++)
          eigvec[j][i]= z[j+i*n];
  }else{
    std::cerr << -info-1 << "[complex] LAPACKE_zheevs: eigval is not converge" << std::endl;
  }
  free(a);
  free(w);
  free(z);
  free(ifail);
}

void eigen_general_solver(const ComplexMatrix &h, ComplexVector &lambda, ComplexMatrix &U_in, ComplexMatrix &V_in, char LOrR='R'){
  /*
   * zgeev computes for an N-by-N complex nonsymmetirc matrix A 
   * the right eigenvertor v(j) of A:
   * A*v(j) = lambda(j) * v(j)
   * The right eigenvector v(j) of A satisfies
   *             A * v(j) = lambda(j) * v(j)
   *  where lambda(j) is its eigenvalue.
   * The left eigenvector u(j) of A satisfies
   *          u(j)**H * A = lambda(j) * u(j)**H
   * where u(j)**H denotes the conjugate transpose of u(j).
   * 
   * V**H is the conjugate transpose of V
   */
  lapack_int n, lda, ldvl, ldvr, info;
  lapack_complex_double *a,*w, *vl,*vr;

  n = h.size();
  lda=n; ldvl=n; ldvr=n;
  
  a = (lapack_complex_double *) calloc(n*n,sizeof(lapack_complex_double));
  w = (lapack_complex_double *) calloc(n,sizeof(lapack_complex_double));
  vl = (lapack_complex_double *) calloc(n*n,sizeof(lapack_complex_double));
  vr = (lapack_complex_double *) calloc(n*n,sizeof(lapack_complex_double));
  
  for(int j=0; j<n; j++)
    for(int i=0; i<n; i++)
      a[i+j*n] = lapack_make_complex_double(h[i][j].real(), h[i][j].imag()); // transform matrix h to col_major

  if(LOrR=='A'){
    info = LAPACKE_zgeev(LAPACK_COL_MAJOR,'N','V',n,a,lda,w,vl,ldvl,vr,ldvr);
    if(info == 0){
      for(int i=0; i<n; i++){
        lambda[i]= w[i];
        for(int j=0; j<n; j++){
          V_in[j][i]= vr[j+i*n];
          U_in[j][i]= vl[j+i*n];
        }
      }
    }else{
      std::cerr << -info-1 << "[Complex] LAPACKE_zgeev: eigval is not converge" << std::endl;
    }
  }else if(LOrR=='L'){
  info = LAPACKE_zgeev(LAPACK_COL_MAJOR,'V','N',n,a,lda,w,vl,ldvl,vr,ldvr);
    if(info == 0){
      for(int i=0; i<n; i++){
        lambda[i]= w[i];
        for(int j=0; j<n; j++)
          U_in[j][i]= vl[j+i*n];
      }
    }else{
      std::cerr << -info-1 << "[Complex] LAPACKE_zgeev: eigval is not converge" << std::endl;
    }
  }else{
    info = LAPACKE_zgeev(LAPACK_COL_MAJOR,'N','V',n,a,lda,w,vl,ldvl,vr,ldvr);
    if(info == 0){
      for(int i=0; i<n; i++){
        lambda[i]= w[i];
        for(int j=0; j<n; j++)
          V_in[j][i]= vr[j+i*n];
      }
    }else{
      std::cerr << -info-1 << "[complex] LAPACKE_zgeev: eigval is not converge" << std::endl;
    }
  }  
  free(a);
  free(w);
  free(vl);
  free(vr);
}

void eigen_general_solver(const DoubleMatrix &h, ComplexVector &lambda, DoubleMatrix &U_in, DoubleMatrix &V_in, char lOrR='R'){
  /*
   * dgeev computes for an N-by-N real nonsymmetirc matrix A 
   * the right eigenvertor v(j) of A:
   * A*v(j) = lambda(j) * v(j)
   * If A is hermite matrix,
   * A=V * lambda * V**H
   * V**H is the conjugate transpose of V
   */
  lapack_int n, lda, ldvl, ldvr, info;
  double *a,*wr, *wi, *vl,*vr;

  n = h.size();
  lda=n; ldvl=n; ldvr=n;
  
  a = (double *) calloc(n*n,sizeof(double));
  wr = (double *) calloc(n,sizeof(double));
  wi = (double *) calloc(n,sizeof(double));
  vl = (double *) calloc(n*n,sizeof(double));
  vr = (double *) calloc(n*n,sizeof(double));
  
  for(int j=0; j<n; j++)
    for(int i=0; i<n; i++)
      a[i+j*n] = h[i][j]; // transform matrix h to col_major

  if(lOrR=='A'){
    info = LAPACKE_dgeev(LAPACK_COL_MAJOR,'V','V',n,a,lda,wr,wi,vl,ldvl,vr,ldvr);
    if(info == 0){
      for(int i=0; i<n; i++){
        lambda[i]= std::complex<double> (wr[i],wi[i]);
        for(int j=0; j<n; j++){
          U_in[j][i]= vl[j+i*n];
          V_in[j][i]= vr[j+i*n];
        }
      }
    }else{
      std::cerr << -info-1 << "[double] LAPACKE_dgeev: eigval is not converge" << std::endl;
    }
  }else if(lOrR=='L'){
    info = LAPACKE_dgeev(LAPACK_COL_MAJOR,'V','N',n,a,lda,wr,wi,vl,ldvl,vr,ldvr);
    if(info == 0){
      for(int i=0; i<n; i++){
        lambda[i]= std::complex<double> (wr[i],wi[i]);
        for(int j=0; j<n; j++)
          U_in[j][i]= vl[j+i*n];
      }
    }else{
      std::cerr << -info-1 << "[double] LAPACKE_dgeev: eigval is not converge" << std::endl;
    }
  }else{
    info = LAPACKE_dgeev(LAPACK_COL_MAJOR,'N','V',n,a,lda,wr,wi,vl,ldvl,vr,ldvr);
    if(info == 0){
      for(int i=0; i<n; i++){
        lambda[i]= std::complex<double> (wr[i],wi[i]);
        for(int j=0; j<n; j++)
          V_in[j][i]= vr[j+i*n];
      }
    }else{
      std::cerr << -info-1 << "[double] LAPACKE_dgeev: eigval is not converge" << std::endl;
    }
  }
  free(a);
  free(wr);
  free(wi);
  free(vl);
  free(vr);
}

void svd_solver(const DoubleMatrix &h_in, DoubleMatrix &u_out, DoubleVector &s_out,DoubleMatrix &vt_out){
  /*
   * dgejsv compute the SVD of a real M-by-N matrix A, where M>=N.
   * A=U*sigma*V^t
   * U is an M-by-M orthonormal matrix
   * V is an N-by-N orthogonal matrix
   * sigma is an M-by-N matrix
   * the diagonal of sigma is stored in sva; they are return in descending order;
   * they are real and non-negative
   *
   * the lapack's dgesvd routine returns V**T, not V
   */
  lapack_int m,n,lda,ldu,ldvt,info;
  double *a,*s,*u,*vt,*superb;
  m=u_out.size(); n=vt_out.size();
  lda=m;ldu=m;ldvt=n;

  a = (double *) calloc(m*n,sizeof(double));  
  s = (double *) calloc(std::min(m,n),sizeof(double));
  u =  (double *) calloc(m*m,sizeof(double));
  vt=  (double *) calloc(n*n,sizeof(double));
  superb = (double *) calloc(std::min(m,n)-1,sizeof(double));

  for(int j=0; j<n; j++)
    for(int i=0; i<m; i++)
      a[i+j*m] = h_in[i][j]; // transform matrix h to col_major

  info=LAPACKE_dgesvd(LAPACK_COL_MAJOR,'A','A',m,n,a,lda,s,u,ldu,vt,ldvt,superb);

  if(info == 0){
    for(int i=0; i<std::min(m,n); i++)
      s_out[i]= s[i];
    for(int i=0; i<m; i++)
      for(int j=0; j<m; j++)
        u_out[j][i]= u[j+i*m];
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        vt_out[j][i]= vt[j+i*n];    
  }else{
    std::cerr << -info-1 << "[double] LAPACKE_dgesvd: svd is false" << std::endl;
  }
  free(a);
  free(s);
  free(u);
  free(vt);
}

void svd_solver(const ComplexMatrix &h_in, ComplexMatrix &u_out, DoubleVector &s_out,ComplexMatrix &vt_out){
  /*
   * dgejsv compute the SVD of a complex M-by-N matrix A, where M>=N.
   * A=U*sigma*V**H
   * U is an M-by-M orthonormal matrix
   * V is an N-by-N orthogonal matrix
   * sigma is an M-by-N matrix
   * the diagonal of sigma is stored in sva; they are return in descending order;
   * they are real and non-negative
   *
   * the lapack's dgesvd routine returns V**H, not V
   */
  lapack_int m,n,lda,ldu,ldvt,info;
  lapack_complex_double *a,*u,*vt;
  double *s, *superb;
  m=h_in.size(); n=h_in[0].size();
  lda=m;ldu=m;ldvt=n;

  a = (lapack_complex_double *) calloc(m*n,sizeof(lapack_complex_double));  
  s = (double *) calloc(std::min(m,n),sizeof(double));
  u =  (lapack_complex_double *) calloc(m*m,sizeof(lapack_complex_double));
  vt=  (lapack_complex_double *) calloc(n*n,sizeof(lapack_complex_double));
  superb = (double *) calloc(std::min(m,n)-1,sizeof(double));

  for(int j=0; j<n; j++)
    for(int i=0; i<m; i++)
      a[i+j*m] = lapack_make_complex_double(h_in[i][j].real(),h_in[i][j].imag());
  // transform matrix h to col_major
  
  info=LAPACKE_zgesvd(LAPACK_COL_MAJOR,'A','A',m,n,a,lda,s,u,ldu,vt,ldvt,superb);

  if(info == 0){
    for(int i=0; i<std::min(m,n); i++)
      s_out[i]= s[i];
    for(int i=0; i<m; i++)
      for(int j=0; j<m; j++)
        u_out[j][i]= u[j+i*m];
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        vt_out[j][i]= vt[j+i*n];
  }else{
    std::cerr << -info-1 << "[complex] LAPACKE_zgesvd: svd is not converge" << std::endl;
  }
  free(a);
  free(s);
  free(u);
  free(vt);
}

/*******************   matrix basic operation  ******************/
void zero_matrix(DoubleMatrix &A){
  for(int i=0; i<A.size(); i++)
    for(int j=0; j<A[0].size(); j++)
      A[i][j]=0;
}

void zero_matrix(ComplexMatrix &A){
  for(int i=0; i<A.size(); i++)
    for(int j=0; j<A[0].size(); j++)
      A[i][j]=0;
}

void one_matrix(DoubleMatrix &A){
  if(A.size()!=A[0].size()){
    std::cerr << "one_matrix: must be square matrix" << std::endl;
    std::exit(1.0);
  }
  for(int i=0; i< A.size(); i++)
    for(int j=0; j<A.size(); j++)
      if(i!=j)
        A[i][j]=0;
      else
        A[i][j]=1.0;
}

void one_matrix(ComplexMatrix &A){
  if(A.size()!=A[0].size()){
    std::cerr << "one_matrix: must be square matrix" << std::endl;
    std::exit(1.0);
  }
  for(int i=0; i< A.size(); i++)
    for(int j=0; j<A.size(); j++)
      if(i!=j)
        A[i][j]=0;
      else
        A[i][j]=1.0;
}

void copy_matrix(DoubleMatrix &A, DoubleMatrix &B){
  if(A.size() != B.size() || A[0].size() != B[0].size() ){
    std::cerr << "copy_matrix: size may be same" << std::endl;
    std::exit(1.0);
  }
  for(int i=0; i<A.size(); i++)
    for(int j=0; j<A[0].size(); j++)
      A[i][j]=B[i][j];
}

void copy_matrix(ComplexMatrix &A, ComplexMatrix &B){
  if(A.size() != B.size() || A[0].size() != B[0].size() ){
    std::cerr << "copy_matrix: size may be same" << std::endl;
    std::exit(1.0);
  }
  for(int i=0; i<A.size(); i++)
    for(int j=0; j<A[0].size(); j++)
      A[i][j]=B[i][j];
}

void copy_diagonal_matrix(DoubleMatrix &A, DoubleVector &B){
  if(A.size() != B.size() || A.size() != A[0].size() ){
    std::cerr << "copy_diagonal_matrix: size may be same" << std::endl;
    std::exit(1.0);
  }
  for(int i=0; i<A.size(); i++)
    for(int j=0; j<A[0].size(); j++)
      if(i!=j)
        A[i][j]=0;
      else
        A[i][j]=B[i];
}

void copy_diagonal_matrix(ComplexMatrix &A, ComplexVector &B){
  if(A.size() != B.size() || A.size() != A[0].size() ){
    std::cerr << "copy_diagonal_matrix: size may be same" << std::endl;
    std::exit(1.0);
  }
  for(int i=0; i<A.size(); i++)
    for(int j=0; j<A[0].size(); j++)
      if(i!=j)
        A[i][j]=0;
      else
        A[i][j]=B[i];
}

void matrix_multiply(DoubleMatrix &A, DoubleMatrix &B, DoubleMatrix &C){
  if(B[0].size() != C.size()){
    std::cerr << "[double] matrix_multiply: col_length of B isn't same as row_length of C" << std::endl;
    exit(1);
  }
  if(A.size() != B.size() || A[0].size() != C[0].size()){
    // std::clog << "[double] matrix_multiply: A's dim doesn't match B and C's" << std::endl;
    A.resize(B.size());
    for(int i=0; i< B.size(); i++)
      A[i].resize(C[0].size());
  }
  for(int i=0; i<B.size(); i++)
    for(int j=0; j<C[0].size(); j++){
      A[i][j]=0;
      for(int k=0; k<C.size(); k++)
        A[i][j] += B[i][k] * C[k][j];
    }
}

void matrix_multiply(ComplexMatrix &A, ComplexMatrix &B, ComplexMatrix &C){
  /*
   * A=B*C
   */
  if(B[0].size() != C.size()){
    std::cerr << "[double] matrix_multiply: col_length of B isn't same as row_length of C" << std::endl;
    exit(1);
  }
  if(A.size() != B.size() || A[0].size() != C[0].size()){
    // std::clog << "[double] matrix_multiply: A's dim doesn't match B and C's" << std::endl;
    A.resize(B.size());
    for(int i=0; i< B.size(); i++)
      A[i].resize(C[0].size());
  }
  for(int i=0; i<B.size(); i++)
    for(int j=0; j<C[0].size(); j++){
      A[i][j]=0;
      for(int k=0; k<C.size(); k++)
        A[i][j] += B[i][k] * C[k][j];
    }
}

void matrix_3multiply(DoubleMatrix &A, DoubleMatrix &B, DoubleMatrix &C, DoubleMatrix &D){
  /*
   * A=B*C*D
   */
  if(B[0].size() != C.size() || C[0].size() != D.size()){
    std::cerr << "[double] matrix_3multiply: col_length of B isn't same as row_length of C" << std::endl;
    exit(1);
  }
  if(A.size() != B.size() || A[0].size() != D[0].size()){
    //std::clog << "[double] matrix_3multiply: A's dim doesn't match B and C's" << std::endl;
    A.resize(B.size());
    for(int i=0; i< B.size(); i++)
      A[i].resize(D[0].size());
  }
  DoubleMatrix temp(B.size());
  for(int i=0; i<B.size(); i++)
    temp[i].resize(C[0].size());
  for(int i=0; i<B.size(); i++)
    for(int j=0; j<C[0].size(); j++){
      temp[i][j]=0;
      for(int k=0; k<C.size(); k++)
        temp[i][j] += B[i][k] * C[k][j];
    }
  for(int i=0; i<B.size(); i++)
    for(int j=0; j<D[0].size(); j++){
      A[i][j]=0;
      for(int k=0; k<D.size(); k++)
        A[i][j] += temp[i][k] * D[k][j];
    }
}

void matrix_3multiply(ComplexMatrix &A, ComplexMatrix &B, ComplexMatrix &C, ComplexMatrix&D){
  /*
   * A=B*C*D
   */
  if(B[0].size() != C.size() || C[0].size() != D.size()){
    std::cerr << "[complex] matrix_3multiply: B, C or D doesn't match" << std::endl;
    exit(1);
  }
  if(A.size() != B.size() || A[0].size() != D[0].size()){
    //std::clog << "[complex] matrix_3multiply: A's dim doesn't match B and C's" << std::endl;
    A.resize(B.size());
    for(int i=0; i< B.size(); i++)
      A[i].resize(D[0].size());
  }
  ComplexMatrix temp(B.size());
  for(int i=0; i<B.size(); i++)
    temp[i].resize(C[0].size());
  for(int i=0; i<B.size(); i++)
    for(int j=0; j<C[0].size(); j++){
      temp[i][j]=0;
      for(int k=0; k<C.size(); k++)
        temp[i][j] += B[i][k] * C[k][j];
    }
  for(int i=0; i<B.size(); i++)
    for(int j=0; j<D[0].size(); j++){
      A[i][j]=0;
      for(int k=0; k<D.size(); k++)
        A[i][j] += temp[i][k] * D[k][j];
    }
}

void matrix_transpose(DoubleMatrix &A){  
  /*
   *A_transpose = A^T
   */
  if(A.size() == A[0].size()){
    double temp;    
    for(int i=0; i< A.size(); i++)
      for(int j=i+1; j < A.size(); j++){
        temp = A[i][j];
        A[i][j] = A[j][i];
        A[j][i] = temp;
      }
  }else{
    DoubleMatrix B(A[0].size());
    for(int i=0; i<B.size(); i++)
      B[i].resize(A.size());
    for(int i=0; i<B.size(); i++)
      for(int j=0; j<B[0].size(); j++)
        B[i][j] = A[j][i];
    A.resize(B.size());
    for(int i=0; i<A.size(); i++)
      A[i].resize(B[0].size());
    for(int i=0; i < B.size(); i++)
      for(int j=0; j < B[0].size(); j++)
        A[i][j] = B[i][j];
  }
}

void matrix_transpose_conjugate(ComplexMatrix &A){  
  /*
   *A_transpose_conjugate = A^\dagger
   */
  if(A.size() == A[0].size()){
    std::complex<double> temp;    
    for(int i=0; i< A.size(); i++){
      A[i][i]=std::conj(A[i][i]);
      for(int j=i+1; j < A.size(); j++){
        temp = A[i][j];
        A[i][j] = std::conj(A[j][i]);
        A[j][i] = std::conj(temp);
      }
    }
  }else{
    ComplexMatrix B(A[0].size());
    for(int i=0; i<B.size(); i++)
      B[i].resize(A.size());
    for(int i=0; i<B.size(); i++)
      for(int j=0; j<B[0].size(); j++)
        B[i][j] = std::conj(A[j][i]);
    A.resize(B.size());
    for(int i=0; i<A.size(); i++)
      A[i].resize(B[0].size());
    for(int i=0; i < B.size(); i++)
      for(int j=0; j < B[0].size(); j++)
        A[i][j] = B[i][j];
  }
}





#endif /* _MATRIX_OPERATION_H_ */
