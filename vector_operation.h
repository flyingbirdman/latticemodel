#ifndef _VECTOR_OPERATION_H_
#define _VECTOR_OPERATION_H_
//#include "nb.h"

void copy_vector(DoubleVector &A, DoubleVector &B){
  /*
   * copy vector B to A.
   */
  if(A.size()!=B.size()){
    std::cerr<<"[double] copy_vector: size is not same" << std::endl;
    std::exit(1.0);
  }
  for(int i=0; i<B.size(); i++)
    A[i]=B[i];
}

void copy_vector(ComplexVector &A, ComplexVector &B){
  /*
   * copy vector B to A.
   */
  if(A.size()!=B.size()){
    std::cerr<<"[complex] copy_vector: size is not same" << std::endl;
    std::exit(1.0);
  }
  for(int i=0; i<B.size(); i++)
    A[i]=B[i];
}

double abs_vector(DoubleVector &A){
  /*
   * ||A||
   */
  double reval=0;
  for(int i=0; i<A.size(); i++)
    reval += A[i]*A[i];
  return sqrt(reval);
}

double abs_vector(ComplexVector &A){
  /*
   * ||A||
   */
  double reval=0;
  for(int i=0; i<A.size(); i++)
    reval += std::norm(A[i]);
  return sqrt(reval);
}

double inner_product(DoubleVector &B, DoubleVector &C){
  /*
   * A= \sum_i B_i*C_i
   */
  if(B.size() != C.size()){
    std::cerr << "dot_product: dim of B and C is not same" << std::endl;
    std::exit(1);
  }

  double A=0;
  for(int i=0; i<B.size(); i++)
    A+=B[i]*C[i];
  return A;
}
std::complex<double> inner_product(ComplexVector &B, ComplexVector &C){
/*
   * A= \sum_i conj(B_i)*C_i
   */
  if(B.size() != C.size()){
    std::cerr << "dot_product: dim of B and C is not same" << std::endl;
    std::exit(1);
  }

  std::complex<double> A(0,0);
  for(int i=0; i<B.size(); i++)
    A+=std::conj(B[i])*C[i];
  return A;
}

void colms_inner_product(DoubleMatrix &A, DoubleMatrix &B, DoubleMatrix &C){
  /*
   * B={B1,B2,...,Bm}
   * C={C1,C2,...,Cn}
   * A= \sum_i B_i*C_i
   */
  int m=B[0].size(); int n=C[0].size();
  if(B.size() != C.size()){
    std::cerr << "inner_product: dim of B and C is not same" << std::endl;
    std::exit(1);
  }
  if(A.size() != m || A[0].size() != n){
    std::cerr << "inner_product: dim of A isn't correct" << std::endl;
    std::exit(1);
  }

  for(int i=0; i<m; i++)
    for(int j=0; j<n; j++){
      A[i][j]=0;
      for(int k=0; k<B.size(); k++)
        A[i][j] += B[k][i]*C[k][j];
    }
}
void colms_inner_product(ComplexMatrix &A, ComplexMatrix &B, ComplexMatrix &C){
  /*
   * B={B1,B2,...,Bm}
   * C={C1,C2,...,Cn}
   * A= \sum_i conj(B_i)*C_i
   */
  int m=B[0].size(); int n=C[0].size();
  if(B.size() != C.size()){
    std::cerr << "inner_product: dim of B and C is not same" << std::endl;
    std::exit(1);
  }
  if(A.size() != m || A[0].size() != n){
    std::cerr << "inner_product: dim of A isn't correct" << std::endl;
    std::exit(1);
  }

  for(int i=0; i<m; i++)
    for(int j=0; j<n; j++){
      A[i][j]=0;
      for(int k=0; k<B.size(); k++)
        A[i][j] += std::conj(B[k][i])*C[k][j];
    }
}
void cross_product(DoubleVector &A, DoubleVector &B, DoubleVector &C){
  /*
   * A= cross(B,C)
   * dim of B and c must be 3
   */
  if(A.size()!=B.size() || B.size() != C.size() || B.size() !=3){
    std::cerr << " cross_product: dim of B and C is not same or dim isn't 3"<<std::endl;
    std::exit(1);
  }
  A[0]=-B[2]*C[1]+B[1]*C[2];
  A[1]=B[2]*C[0]-B[0]*C[2];
  A[2]=-B[1]*C[0]+B[0]*C[1];
}
double triple_product(DoubleVector &A, DoubleVector &B, DoubleVector &C){
  /*
   * dot(A, cross(B,C))
   * all of dims must be 3
   */
  if(A.size()!=B.size() || B.size()!=C.size() || B.size()!=3){
    std::cerr << " cross_product: dim of A, B and C is not same or dim isn't 3"<<std::endl;
    std::exit(1);
  }
  return A[0]*(-B[2]*C[1]+B[1]*C[2]) + A[1]*(B[2]*C[0]-B[0]*C[2]) + A[2]*(-B[1]*C[0]+B[0]*C[1]);
}

double norm2_vector(DoubleVector &A){
  /*
   * sqrt{sum_i |A_i|^2 }
   */
  double norm2=0;
  for(int i=0; i<A.size(); i++)
    norm2 += A[i]*A[i];
  return sqrt(norm2);
}
double norm2_vector(ComplexVector &A){
  /*
   * sqrt{sum_i |A_i|^2 }
   */
  double norm2=0;
  for(int i=0; i<A.size(); i++)
    norm2 += A[i].real()*A[i].real()+A[i].imag()*A[i].imag();
  return sqrt(norm2);
}

/******************   matrix and vector   ******************/
template <class T>
void out_product(std::vector<std::vector<T> > &H, std::vector<T> &V1, std::vector<T> &V2){
  if(H.size() != V1.size() || H[0].size() != V2.size()){
    H.resize(V1.size());
    for(int i=0; i<V1.size(); i++)
      H[i].resize(V2.size());
  }
  for(int i=0; i<V1.size(); i++)
    for(int j=0; j<V2.size(); j++)
      H[i][j]=V1[i]*V2[j];
}

template <class T1, class T2, class T3>
void matrix_vector_product(std::vector<T1> &A, std::vector<std::vector<T2> > &B, std::vector<T3> &C){
  /*
   * [Vector]A = [Matrix]B * [Vector]C
   */
  if(B[0].size() != C.size()){
    std::cerr << "[double]matrix_vector_multiply: dim of B and C not match" << std::endl;
    std::exit(1);
  }
  for(int i=0; i<B.size(); i++){
    A[i]=0;
    for(int j=0; j<B[0].size(); j++){
      A[i] += B[i][j]*C[j];
    }
  }
}

double inner_product_with_matrix(DoubleVector &V1, DoubleMatrix &H, DoubleVector &V2){  
  /*
   * calculate: <V1| H |V2>
   */
  if(V1.size() != H.size() || H[0].size() != V2.size()){
    std::cerr << "[double] coupling_or_expectation: dim of V1, H and V2 don't match" << std::endl;
    std::exit(1);
  }
  DoubleVector temp_vec(H.size());
  for(int i=0; i< H.size(); i++){
    temp_vec[i]=0;
    for(int j=0; j< H[0].size(); j++)
      temp_vec[i] += H[i][j]*V2[j];
  }

  double retval=0;
  for(int i=0; i<H.size(); i++)
    retval += V1[i]*temp_vec[i];
  return retval;
}

std::complex<double> inner_product_with_matrix(ComplexVector &V1, ComplexMatrix &H, ComplexVector &V2){
  /*
   * calculate: <V1| H |V2>
   */
  if(V1.size() != H.size() || H[0].size() != V2.size()){
    std::cerr << "[double] coupling_or_expectation: dim of V1, H and V2 don't match" << std::endl;
    std::exit(1);
  }
  ComplexVector temp_vec(H.size());
  for(int i=0; i< H.size(); i++){
    temp_vec[i]=0;
    for(int j=0; j< H[0].size(); j++)
      temp_vec[i] += H[i][j]*V2[j];
  }
  std::complex<double> retval(0,0);
  for(int i=0; i<H.size(); i++)
    retval += std::conj(V1[i])*temp_vec[i];
  return retval;
}


void colms_inner_product_with_matrix(DoubleMatrix &V12, DoubleMatrix &V1, DoubleMatrix &H, DoubleMatrix &V2){  
  /*
   * v1={v11,v12,v13...v1m}
   * v2={v21,v22,v23...v2n}
   * output is m*n matrix with element <v1i|H|v2j>
   * calculate: <V1| H |V2>
   */
  int m=V1[0].size(), n=V2[0].size();
  if(V1.size() != H.size() || H[0].size() != V2.size()){
    std::cerr << "[double] colms_inner_product_matrix: dim of V1, H and V2 don't match" << std::endl;
    std::exit(1);
  }
  if(V12.size() != m || V12[0].size() != n){
    std::cerr << "[double] colms_inner_product_matrix: dim of V12 isn't correct" << std::endl;
    std::exit(1);
  }
  DoubleMatrix temp_vec(H.size());
  for(int i=0; i<H.size(); i++)
    temp_vec[i].resize(n);
  for(int i=0; i< H.size(); i++)
    for(int j=0; j< n; j++){
      temp_vec[i][j]=0;
      for(int k=0; k< H[0].size(); k++)
        temp_vec[i][j] += H[i][k]*V2[k][j];  
    }

  for(int i=0; i<m; i++)
    for(int j=0; j<n;j++){
      V12[i][j]=0;
      for(int k=0; k<H.size(); k++)
        V12[i][j] += V1[k][i]*temp_vec[k][j];
    }
}

void colms_inner_product_with_matrix(ComplexMatrix &V12, ComplexMatrix &V1, ComplexMatrix &H, ComplexMatrix &V2){
  /*
   * v1={v11,v12,v13...v1m}
   * v2={v21,v22,v23...v2n}
   * output is m*n matrix with element <v1i|H|v2j>
   * calculate: <V1| H |V2>
   */
  int m=V1[0].size(), n=V2[0].size();
  if(V1.size() != H.size() || H[0].size() != V2.size()){
    std::cerr << "[double] colms_inner_product_matrix: dim of V1, H and V2 don't match" << std::endl;
    std::exit(1);
  }
  if(V12.size() != m || V12[0].size() != n){
    std::cerr << "[double] colms_inner_product_matrix: dim of V12 isn't correct" << std::endl;
    std::exit(1);
  }
  ComplexMatrix temp_vec(H.size());
  for(int i=0; i<H.size(); i++)
    temp_vec[i].resize(n);
  for(int i=0; i< H.size(); i++)
    for(int j=0; j< n; j++){
      temp_vec[i][j]=0;
      for(int k=0; k< H[0].size(); k++)
        temp_vec[i][j] += H[i][k]*V2[k][j];  
    }

  for(int i=0; i<m; i++)
    for(int j=0; j<n;j++){
      V12[i][j]=0;
      for(int k=0; k<H.size(); k++)
        V12[i][j] += std::conj(V1[k][i])*temp_vec[k][j];
    }
}

template <class T>
void get_colm_matrix_i(std::vector<T> &V, std::vector<std::vector<T> > &H, int k){
  /*
   * get the column of matrix H, V = H[:,k];
   */
  if(V.size() != H.size())
    V.resize(H.size());
  for(int i=0; i<H.size(); i++)
    V[i]=H[i][k];
}

#endif /* _VECTOR_OPERATION_H_ */
    
