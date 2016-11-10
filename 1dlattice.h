#ifndef _1DLATTICE_H_
#define _1DLATTICE_H_

#include "nb.h"

/******************  1D Mesh and DOS  ******************/
struct D1MeshPt{
  DoubleVector en; // energy
  DoubleVector d1; // energy derivative along kx
  // DoubleMatrix ev; // eigenvector of energy, en
};
template<class T>
class D1Mesh{
  /*
   * calculate En(kx,ky)
   * fourier transform of En(kx,ky)
   * 
   */
 public:
  double px1, px2;
  int Lq, numOfBand;
  double eps;
  Eigen<T> *esol;  
  std::vector<D1MeshPt > grid;
  D1Mesh(Eigen<T> *eig, double x1, double x2, int b, int L=128,double dk=0.00001):
      esol(eig), px1(x1),px2(x2),numOfBand(b),Lq(L),eps(dk){
    // currently, Lq must be equal to Lq;
    pt.en.resize(numOfBand);
    pt.d1.resize(numOfBand);
    ki.resize(1);
    veci.resize(esol->hamiltonian.size());
    for(int i=0; i<veci.size(); i++)
      veci[i].resize(numOfBand);
  }
  D1Mesh(std::string &infile){
    std::ifstream in(infile.c_str());
    if(!in){
    std::cerr << "[D1Mesh] error: can't open config file" << std::endl;
    std::exit(1.0);
    }
    datFromFile(in);
  }
  
  double kx_i(int index);
  void cal_grid();
  void d1FourTransform(int i=0, std::ostream &out=std::cout);  // only fourier a single band
  void datToFile(std::ostream &out);
  void datFromFile(std::istream &in);
  
  DoubleVector ki;
  std::vector<std::vector<T> > veci;
  D1MeshPt pt;
};

template<class T>
class D1Dos{
  /*
   * DOS 
   *
   */   
 public:
  D1Mesh<T>* d1m;
  DoubleMatrix xx;
  DoubleMatrix yy;
  double emin, emax, zoneArea;
  int ne, nq, numOfBand;
  double de, dq;
  double qx1, qx2;
  double eps; 
  D1Dos(D1Mesh<T>* dm, double e1=0.0, double e2=8.0, double n=1000,double tol=0.0000001):
      d1m(dm),zoneArea(dm->px2-dm->px1),emin(e1),emax(e2),ne(n),eps(tol){
    nq=d1m->Lq;
    qx1=d1m->px1; qx2=d1m->px2;
    dq=(qx2-qx1)/nq;
    de=(emax-emin)/ne;
    numOfBand=d1m->numOfBand;
    xx.resize(ne);
    yy.resize(ne);
    for(int i=0; i<ne; i++){
      xx[i].resize(1+numOfBand);
      yy[i].resize(1+numOfBand);
    }
  }
  void dosFunc();
  void hilbertTransform();
  void xyToFile(std::ostream &out);
  // double (*MatrixElement)(DoubleVector &k1, DoubleVector &k2);
  // only for correlation function
};

#endif /* _1DLATTICE_H_ */
