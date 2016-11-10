#ifndef _NRLIB_H_
#define _NRLIB_H_

#include "include/nr3.h"
#include "include/qgaus.h"
#include "include/quad3d.h"
#include "include/interp_1d.h"
#include "include/quadrature.h"
#include "include/romberg.h"
#include "include/fourier.h"
#include "include/fourier_ndim.h"

#ifndef _VECTOR_MATRIX_DYPEDEF_
#define _VECTOR_MATRIX_DYPEDEF_
typedef std::vector<double> DoubleVector;
typedef std::vector<std::complex<double> > ComplexVector;
typedef std::vector<std::vector<double> > DoubleMatrix;
typedef std::vector<std::vector<std::complex<double> > > ComplexMatrix;
typedef std::vector<std::vector<std::vector<double> > > DoubleTensor;
typedef std::vector<std::vector<std::vector<std::complex<double> > > > ComplexTensor;
#endif /* _VECTOR_MATRIX_DYPEDEF_*/

/******************  fourier transformation  ******************/
void fft_1d(DoubleVector &h_k, int isign)
{
  /*
   * H_n = \sum_{k=0}^{N-1} h_k exp{isign * 2\pi ikn/N}
   * isign: the sign before 'i', +1 or -1
   * [p2x]: isign= 1;
   * [t2e]: isign= 1;
   */
  VecDoub lh_k(h_k.size());
  for(int i=0; i<lh_k.size(); i++)
    lh_k[i]=h_k[i];
  four1(lh_k, isign);
  for(int i=0; i < h_k.size(); i++)
    h_k[i] = h_k[i]/(h_k.size()/2.0-1.0);
}

void fft_1d(ComplexVector &h_k, int isign)
{
  /*
   * H_n = \sum_{k=0}^{N-1} h_k exp{isign * 2\pi ikn/N}
   * isign: the sign before 'i', +1 or -1
   * [p2x]: isign= 1;
   * [t2e]: isign= 1;
   */  
  VecDoub lh_k(h_k.size()*2);
  for(int i=0; i<h_k.size(); i++){
    lh_k[2*i]=h_k[i].real();
    lh_k[2*i+1]=h_k[i].imag();
  }
  four1(lh_k, isign);
  double Delta = h_k.size()-1.0;
  for(int i=0; i<h_k.size(); i++)
    h_k[i]=std::complex<double> (lh_k[2*i]/Delta,lh_k[2*i+1]/Delta);
}

void fft_2d(ComplexMatrix &h_k, int isign){
  /*
   * H_n = \sum_{k1,k2=0}^{N1-1,N2-1} h_k1_k2
   *           * exp{isign * 2\pi ik1 n/N1}
   *           * exp{isign * 2\pi ik2 n/N2};
   * isign: the sign before 'i', +1 or -1
   * [p2x]: isign= 1;
   * [t2e]: isign= 1;
   * data: 1d-array data[0 ... (2N1N2-1)]
   */
  int m=h_k.size();
  int n=h_k[0].size();
  Doub *lh_k;
  lh_k = (Doub*) calloc((2*m*n-1),sizeof(Doub));
  for(int i=0; i<m; i++)
    for(int j=0; j<n; j++){
      lh_k[2*(j+i*n)]=h_k[i][j].real();
      lh_k[2*(j+i*n)+1]=h_k[i][j].imag();
    }
  VecInt nn(2);
  nn[0]=m; nn[1]=n;
  fourn(lh_k,nn,isign);
  double Delta= (m-1.0)*(n-1.0);
  for(int i=0; i<m; i++)
    for(int j=0; j<n; j++)
      h_k[i][j]=std::complex<double> (lh_k[2*(j+i*n)]/Delta, lh_k[2*(j+i*n)+1]/Delta);
  free(lh_k);
}

void print_fft_1d(ComplexVector &h_k, double delta, std::ostream &out=std::cout){
  out << "# fft_1d" << std::endl;
  double n=h_k.size();
  for(int i=1;i<n/2;i++)
    out<<-1.0/n*(n/2-i)/delta << "  " << abs(h_k[n/2+i]) << "  " << imag(h_k[n/2+i]) << std::endl;
  for(int i=0;i<n/2;i++)
    out<<1.0/n*i/delta << "  " << abs(h_k[i]) << "  " << imag(h_k[i]) << std::endl;
  out << std::endl;
}

void print_fft_2d(ComplexMatrix &h_k, double delta1, double delta2, std::ostream &out=std::cout){
  /*
   *delta: time-interval of the discrete data in the [t-E space]
   *       or space-interval of discrete data in the [r-p space]
   */
  out << "# fft_2d" << std::endl;
  double n1=h_k.size(), n2= h_k[0].size();
  for(int i=1; i< n1/2; i++){
    for(int j=1; j< n2/2; j++){
      out<< -1.0/n1*(n1/2-i)/delta1 << "  " << -1.0/n2*(n2/2-j)/delta2 << "  " << abs(h_k[n1/2+i][n2/2+j]) << "  " << imag(h_k[n1/2+i][n2/2+j]) <<std::endl;
    }
    for(int j=0; j< n2/2; j++){
      out<< -1.0/n1*(n1/2-i)/delta1 << "  " << 1.0/n2*j/delta2 << "  " << abs(h_k[n1/2+i][j]) << "  " << imag(h_k[n1/2+i][j]) <<std::endl;
    }
    out << std::endl;
  }
  for(int i=0; i< n1/2; i++){
    for(int j=1; j< n2/2; j++){
      out<< 1.0/n1*i/delta1 << "  " << -1.0/n2*(n2/2-j)/delta2 << "  " << abs(h_k[i][n2-j]) << "  " << imag(h_k[i][n2/2+j]) <<std::endl;
    }
    for(int j=0; j< n2/2; j++){
      out<< 1.0/n1*i/delta1 << "  " << 1.0/n2*j/delta2 << "  " << abs(h_k[i][j]) << "  " << imag(h_k[i][j]) <<std::endl;
    }
    out << std::endl;
  }
  out << std::endl;
}
void fft_3d(ComplexTensor &h_k, int isign){
  /*
   * H_n = \sum_{k1,k2,k3=0}^{N1-1,N2-1,N3-1} h_k1_k2_k3
   *           * exp{isign * 2\pi ik1 n/N1}
   *           * exp{isign * 2\pi ik2 n/N2}
   *           * exp{isign * 2\pi ik3 n/N3};
   * isign: the sign before 'i', +1 or -1
   * [p2x]: isign= 1;
   * [t2e]: isign= 1;
   * data: 1d-array data[0 ... (2N1N2N3-1)]
   */
      
  int n1=h_k.size();
  int n2=h_k[0].size();
  int n3=h_k[0][0].size();
  Doub* lh_k;
  lh_k = (Doub*) calloc((2*n1*n2*n3-1),sizeof(Doub));
  for(int i=0; i<n1; i++)
    for(int j=0; j<n2; j++)
      for(int k=0; k<n3; k++){
        lh_k[2*(k+j*n3+i*n3*n2)]=h_k[i][j][k].real();
        lh_k[2*(k+j*n3+i*n3*n2)+1]=h_k[i][j][k].imag();
      }
  VecInt nn(3);
  nn[0]=n1; nn[1]=n2; nn[2]=n3;
  fourn(lh_k,nn,isign);
  double Delta= (n1-1.0)*(n2-1.0)*(n3-1.0);
  for(int i=0; i<n1; i++)
    for(int j=0; j<n2; j++)
      for(int k=0; k<n3; k++)
        h_k[i][j][k]=std::complex<double> (lh_k[2*(k+j*n3+i*n3*n2)]/Delta, lh_k[2*(k+j*n3+i*n3*n2)+1]/Delta);
  free(lh_k);
}

/******************  2d quadrature  ******************/
struct NRf2_2d{
  Doub xsav;
  Doub (*func2d)(const Doub, const Doub);
  Doub operator()(const Doub y){
    return func2d(xsav,y);
  }
};
struct NRf1_2d{
  NRf2_2d f2;
  Doub (*y1)(Doub);
  Doub (*y2)(Doub);
  NRf1_2d(Doub yy1(Doub), Doub yy2(Doub)): y1(yy1), y2(yy2){};  // 
  Doub operator()(const Doub x){
    f2.xsav=x;
    return qgaus(f2,y1(x),y2(x));
  } 
};
template <class T>
Doub quad2d(T &func, const Doub x1, const Doub x2, Doub y1(Doub), Doub y2(Doub)){
  /* 
   * minium x1 and maximum x2, low_boundary y1(x1) and up_boundary y2(x)
   * must be given.
   */
  NRf1_2d f1(y1,y2);
  f1.f2.func2d = func;
  return qgaus(f1,x1,x2);
}

#endif /* _NRLIB_H_ */
