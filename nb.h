/**
 * @author: Kai Yang
 * @mail: yangkai0208@gmail.com
 *
 * @brief calculate the general property of band structure
 * in optical lattice, i.e, DOS, energy spectrum, maxmium localized
 * wannier function. Based on the basic property, the multi-band
 * hamiltonian in the quasi-k space and topological property for some 
 * isolated or composite bands is also analysable.
 * 
 * @date: June 7, 2015
 * 
 **/
#ifndef _NB_H_
#define _NB_H_
// numerical calculation of band structure

// all the #include's for numeric calculation
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iterator>
#include <vector>
#include <set>
#include <list>
#include <utility> //pair struct
#include <algorithm>
#define pi 3.141592653589793
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
lapack_complex_float lapack_make_complex_float( float re, float im );
lapack_complex_double lapack_make_complex_double( double re, double im );
#include <lapacke/lapacke.h>

lapack_complex_float lapack_make_complex_float( float re, float im ){
  return std::complex<float> (re,im);
}
lapack_complex_double lapack_make_complex_double( double re, double im ){
  return std::complex<double> (re,im);
}

#ifndef _VECTOR_MATRIX_DYPEDEF_
#define _VECTOR_MATRIX_DYPEDEF_
typedef std::vector<double> DoubleVector;
typedef std::vector<std::complex<double> > ComplexVector;
typedef std::vector<std::vector<double> > DoubleMatrix;
typedef std::vector<std::vector<std::complex<double> > > ComplexMatrix;
typedef std::vector<std::vector<std::vector<double> > > DoubleTensor;
typedef std::vector<std::vector<std::vector<std::complex<double> > > > ComplexTensor;
#endif /* _VECTOR_MATRIX_DYPEDEF_*/

/*******************    adjacent point    ******************/
template<class T=double>
struct Adj_list{
  DoubleVector vec;
  std::list<std::pair<int,T> > indx;//k_index and coupling strength
  Adj_list(DoubleVector &tmp_vec, std::list<std::pair<int,double> > &tmp_list): vec(tmp_vec), indx(tmp_list){
  };
};
typedef std::vector<Adj_list<double> > Adj_listVector;
typedef std::vector<Adj_list<std::complex<double> > > cAdj_listVector;
/*******************   vector compare criterion   ******************/
class VecCmp{ //define the compare criterion for k points 
 public:
  enum cmp_mode {normal, reverse};
  VecCmp(cmp_mode m=VecCmp::normal) : mode(m){
  }
  VecCmp(const VecCmp &cmp){
    mode=cmp.mode;
  }
  bool operator() (const DoubleVector &v1, const DoubleVector &v2) const;
  bool operator== (const VecCmp& vc){
    return mode == vc.mode;
  }
 private:
  cmp_mode mode;
};
typedef std::set<DoubleVector, VecCmp> cmpSet;  

bool VecCmp::operator() (const DoubleVector &v1, const DoubleVector &v2) const {
  double nv1=0; double nv2=0; bool retbool;
  for(int i=0; i < v1.size(); i++){
    nv1 += v1[i]*v1[i];
    nv2 += v2[i]*v2[i];
  }
  if(fabs(nv1-nv2) > 10e-14)
    return (mode == VecCmp::normal ? nv1 < nv2 : nv1 > nv2);
  for(int i=0; i < v1.size(); i++){
    if(fabs(v1[i]-v2[i])>10e-14)
      return (mode == VecCmp::normal ? v1[i] < v2[i] : v1[i] > v2[i]);
    else
      retbool=false;
  }
  return retbool;
}

/*******************   lattice struct, adj_k & adj_k_network   ******************/
struct LatticeStruct{
  DoubleMatrix basis;    // basis[i][:]
  DoubleVector coupling; // link to basis[i]
  int dim;               // dim of real space
  int numOfBasis;
  int numOfMagnetic;     // total # of quantum magnetic index
  std::string struct_name;
  std::string socOrNot;
  DoubleMatrix Sx, Sz;
  ComplexMatrix Sy;
  double Qx, Qy, Qz, Omega;  // (k-QxSx)**2+(ky-QySy)**2+(kz-QzSz)**2+Omega*Sz
  void generate_spin_matrix();
  LatticeStruct(DoubleMatrix &b, DoubleVector &c, int d, std::string &name,
                std::string soc="no", int tm=1, double omg=0, double qx=0, double qy=0, double qz=0):
      basis(b),coupling(c),struct_name(name),dim(d),numOfBasis(b.size()),
      socOrNot(soc),numOfMagnetic(tm),Omega(omg),Qx(qx),Qy(qy),Qz(qz){
    if(dim!=basis[0].size()){
      std::cerr << "LatticeStruct: wrong dimension !" << std::endl;
      std::exit(1.0);
    }
    if(socOrNot=="yes" || socOrNot=="y"){
      if(numOfMagnetic>1)
        generate_spin_matrix();
      else{
        std::cerr << "[LatticeStruct]: soc==yes, numOfMagnetic>1, pls check config!" << std::endl;
        std::exit(1.0);
      }
    }
  }
};
struct Adj_k{
  LatticeStruct ls;
  cmpSet k_set;
  double k_max;
  Adj_k(LatticeStruct &s, double max=10):ls(s),k_max(max){};
  int dimOfHilbert;
  double dist;  // just a temp variable;
  double length(DoubleVector &v){
    dist=0;
    for(int i=0; i<v.size(); i++)
      dist+=v[i]*v[i];
    dist=sqrt(dist);
    return dist;
  }
  bool inRange(DoubleVector &v){  
    return length(v) < k_max;
  }
  void print_k_set();
  void generate_adj_k();  
};

struct Adj_k_network{
  // define a network of k space moentum transfer
  Adj_k_network(LatticeStruct &s):ls(s),dim(s.dim){ };
  LatticeStruct ls;
  Adj_listVector adj_network;
  void get_adj_network(Adj_k &ak);
  void print_adj_network();
  int dimOfHilbert;
  int dim;  // dim of real space
};
/*******************   eigen system for band structure  ******************/
void get_off_diag_hamiltonian(DoubleMatrix &h, const  Adj_k_network &akn);
void get_off_diag_hamiltonian(ComplexMatrix &h, const  Adj_k_network &akn);
void get_diag_hamiltonian(const DoubleVector &p, DoubleMatrix &h, const Adj_k_network &akn);
void get_diag_hamiltonian(const DoubleVector &p, ComplexMatrix &h, const Adj_k_network &akn);

template<class T>
class Eigen{
 public:
  Eigen(Adj_k_network *nw, double fft=128, int n=8):
      akn(nw),
      dim(nw->dim),
      fftsize(fft),
      numOfBand(n),
      eigval(DoubleVector(n)){
    if(typeid(T).name()==typeid(double).name())
      data_type="real";
    else if(typeid(T).name()==typeid(std::complex<double>).name())
      data_type="complex";
    else{
      std::cerr<< "Eigen: unknown data type" << std::endl;
      std::exit(1.0);
    }
    hamiltonian.resize(akn->dimOfHilbert);
    eigvec.resize(akn->dimOfHilbert);
    for(int i=0; i<nw->dimOfHilbert; i++)
      hamiltonian[i].resize(akn->dimOfHilbert);
    for(int i=0; i<nw->dimOfHilbert; i++)
      eigvec[i].resize(n);
    get_off_diag_hamiltonian(hamiltonian,*akn);
  };
  std::string data_type;
  std::vector<std::vector<T> > hamiltonian;
  Adj_k_network *akn;  
  void cal_ei(DoubleVector &ki, DoubleVector &val, std::vector<std::vector<T> > &vec, char jobz='N');
  // api for other computation
  void cal_loop(DoubleMatrix &loop, std::ofstream &out);
  void cal_d1zone(DoubleMatrix &zone, std::ofstream &out);
  void cal_d2zone(DoubleMatrix &zone, std::ofstream &out);
  void cal_d3zone(DoubleMatrix &zone, std::ofstream &out);
  void obtain_eigval_dat(const DoubleVector &k, std::ofstream &out);  
  
  int dim;
  int numOfBand;
  int fftsize;
  DoubleVector eigval;
  std::vector<std::vector<T> > eigvec;
};

/*******************   subspace analysis    ******************/
template <class T>
class SubEigen{
 public:
  std::string data_type;
  std::vector<std::vector<T> > h1, h2, lh;
  Adj_k_network *akn1, *akn2;
  int bl, bu;  // lowest and highest band: define the subspace
  
  void obtain_subspace_hamiltonian(DoubleVector &ki, std::ofstream &out);
  void obtain_subspace_eigenvector(DoubleVector &ki, std::ostream &out);
  void order_v1_basis();
  void gauge_fix();
  void cal_lk_subspace(DoubleVector &ki);
  void cal_d1_subspace(DoubleMatrix &zone, std::ofstream &out);
  void cal_d2_subspace(DoubleMatrix &zone, std::ofstream &out);
  void cal_d3_subspace(DoubleMatrix &zone, std::ofstream &out);

  std::vector<int> v1_order;
  int dim;
  int fftsize;
  int dimOfHilbert;
  int dimOfSub;
  DoubleVector val1, val2;
  std::vector<std::vector<T> > vec1, vec2;  // vec1 is the basis
  std::vector<std::vector<T> > lvec;

  SubEigen(std::string t, Adj_k_network *nwa, Adj_k_network *nwb, int b1, int b2,int fft=128):
      akn1(nwa), akn2(nwb), bl(b1), bu(b2), fftsize(fft),
      dim(nwa->dim), dimOfHilbert(nwa->dimOfHilbert),dimOfSub(b2-b1+1){
    if(typeid(T).name()==typeid(double).name())
      data_type="real";
    else if(typeid(T).name()==typeid(std::complex<double>).name())
      data_type="complex";
    else{
      std::cerr<< "SubEigen: unknown data type" << std::endl;
      std::exit(1.0);
    }
    if(akn1->dimOfHilbert != akn2->dimOfHilbert){
      std::cerr << "SubEigen: dim1 isn't same as dim2" << std::endl;
      std::exit(1.0);
    }
    if(bl>bu){
      std::cerr << "SubEigen: dim of subspace less than 0" << std::endl;
      std::exit(1.0);
    }
    h1.resize(dimOfHilbert);
    h2.resize(dimOfHilbert);
    for(int i=0; i<dimOfHilbert; i++){
      h1[i].resize(dimOfHilbert);
      h2[i].resize(dimOfHilbert);
    }
    vec1.resize(dimOfHilbert);
    vec2.resize(dimOfHilbert);
    get_off_diag_hamiltonian(h1,*akn1);
    get_off_diag_hamiltonian(h2,*akn2);  
    for(int i=0; i< dimOfHilbert; i++){
      vec1[i].resize=bu-bl+1;
      vec1[i].resize=bu-bl+1;
    }
    val1.resize(bu-bl+1);
    val2.resize(bu-bl+1);
    lh.resize(bu-bl+1);
    lvec.resize(bu-bl+1);
    for(int i=0; i<bu-bl+1; i++){
      lh[i].resize(bu-bl+1);
      lvec[i].resize(bu-bl+1);
    }
    v1_order.resize(bu-bl+1);
  }
};

/*******************   print matrix and vector  ******************/
void print_matrix(DoubleMatrix &A, std::ostream &out=std::clog){
  for(int i=0; i<A.size(); i++){
    for(int j=0; j<A[i].size(); j++)
      out << A[i][j] << "     ";
    out << std::endl;
  }
  out << std::endl;
}
void print_matrix(ComplexMatrix &A, std::ostream &out=std::clog){
  for(int i=0; i<A.size(); i++){
    for(int j=0; j<A[i].size(); j++)
      out << A[i][j] << "     ";
    out << std::endl;
  }
  out << std::endl;
}
void print_vector(DoubleVector &V, std::ostream &out=std::clog){
  for(int i=0; i<V.size(); i++)
    out << V[i] << "     ";
  out << std::endl;
  out << std::endl;
}
void print_vector(ComplexVector &V, std::ostream &out=std::clog){
  for(int i=0; i<V.size(); i++)
    out << V[i] << "     ";
  out << std::endl;
  out << std::endl;
}

void print_matrix_m(ComplexMatrix &A, std::ostream &out=std::clog){
  out.setf(std::ios::fixed);
  out.precision(5);
  out << "{";
  for(int i=0; i<A.size(); i++){
    out << "{";
    for(int j=0; j<A[i].size(); j++){
      if(j!=A[i].size()-1)
        out << A[i][j].real() << "+" << A[i][j].imag() <<"I,   ";
      else
        out << A[i][j].real() << "+" << A[i][j].imag() << "I ";
    }
    if(i!=A.size()-1)
      out << "}," << std::endl;
    else
      out <<"}" << std::endl;
  }
  out << "};" << std::endl;  // 
}
void print_matrix_m(DoubleMatrix &A, std::ostream &out=std::clog){  
  out.setf(std::ios::fixed);
  out.precision(5);
  out << "{";
  for(int i=0; i<A.size(); i++){
    out << "{";
    for(int j=0; j<A[i].size(); j++){
      if(j!=A[i].size()-1)
        out << A[i][j] <<",  ";
      else
        out << A[i][j];
    }
    if(i!=A.size()-1)
      out << "}," << std::endl;
    else
      out <<"}" << std::endl;
  }
  out << "};" << std::endl;
}
void print_vector_m(DoubleVector &V, std::ostream &out=std::clog){
  out.setf(std::ios::fixed);
  out.precision(5);
  out << "{";
  for(int i=0; i<V.size(); i++){
    if(i!=V.size()-1)
      out << V[i] << ",    ";
    else
      out << V[i] << " ";
  }
  out << "};" << std::endl;
  out << std::endl;
}
void print_vector_m(ComplexVector &V, std::ostream &out=std::clog){
  out.setf(std::ios::fixed);
  out.precision(5);
  out << "{";
  for(int i=0; i<V.size(); i++){
    if(i!=V.size()-1)
      out << V[i].real() <<"+"<< V[i].imag() << "I,   ";
    else
      out << V[i].real() <<"+"<< V[i].imag() << "I ";
  }
  out << "};" << std::endl;
  out << std::endl;
}

void print_matrix_colm_order_m(ComplexMatrix &A, std::vector<int> &ord, std::ostream &out=std::clog){
  if(ord.size()!=A[0].size()){
    std::cerr << "print_matrix_colm_oder_m: size of order problem. " << std::endl;
    std::exit(1.0);
  }
  out.setf(std::ios::fixed);
  out.precision(5);
  out << "{";
  for(int i=0; i<A.size(); i++){
    out << "{";
    for(int j=0; j<A[i].size(); j++){
      if(j!=A[i].size()-1)
        out << A[i][ord[j]].real() << "+" << A[i][ord[j]].imag() <<"I,   ";
      else
        out << A[i][ord[j]].real() << "+" << A[i][ord[j]].imag() << "I ";
    }
    if(i!=A.size()-1)
      out << "}," << std::endl;
    else
      out <<"}" << std::endl;
  }
  out << "};" << std::endl;  // 
}
void print_matrix_colm_order_m(DoubleMatrix &A, std::vector<int> &ord, std::ostream &out=std::clog){
  if(ord.size()!=A[0].size()){
    std::cerr << "print_matrix_colm_oder_m: size of order problem. " << std::endl;
    std::exit(1.0);
  }
  out.setf(std::ios::fixed);
  out.precision(5);
  out << "{";
  for(int i=0; i<A.size(); i++){
    out << "{";
    for(int j=0; j<A[i].size(); j++){
      if(j!=A[i].size()-1)
        out << A[i][ord[j]] <<",  ";
      else
        out << A[i][ord[j]];
    }
    if(i!=A.size()-1)
      out << "}," << std::endl;
    else
      out <<"}" << std::endl;
  }
  out << "};" << std::endl;
}

void print_vector_order_m(DoubleVector &V, std::vector<int> &ord, std::ostream &out=std::clog){
  if(ord.size() != V.size()){
    std::cerr << "print_vector_order_m: size of order problem. " << std::endl;
    std::exit(1.0);
  }
  out.setf(std::ios::fixed);
  out.precision(5);
  out << "{";
  for(int i=0; i<V.size(); i++){
    if(i!=V.size()-1)
      out << V[ord[i]] << ",    ";
    else
      out << V[ord[i]] << " ";
  }
  out << "};" << std::endl;
  out << std::endl;
}
void print_vector_order_m(ComplexVector &V, std::vector<int> &ord, std::ostream &out=std::clog){
  if(ord.size() != V.size()){
    std::cerr << "print_vector_order_m: size of order problem. " << std::endl;
    std::exit(1.0);
  }
  out.setf(std::ios::fixed);
  out.precision(5);
  out << "{";
  for(int i=0; i<V.size(); i++){
    if(i!=V.size()-1)
      out << V[ord[i]].real() <<"+"<< V[ord[i]].imag() << "I,   ";
    else
      out << V[ord[i]].real() <<"+"<< V[ord[i]].imag() << "I ";
  }
  out << "};" << std::endl;
  out << std::endl;
}

#endif /* _NB_H */
