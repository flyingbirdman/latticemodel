#ifndef _2DLATTICE_H_
#define _2DLATTICE_H_

#include "nb.h"

/******************  2D Mesh and DOS  ******************/
struct D2MeshPt{
  DoubleVector en; // energy
  DoubleVector d1; // energy derivative along kx
  DoubleVector d2; // energy derivative along ky
  // DoubleMatrix ev; // eigenvector of energy, en
};
template<class T>
class D2Mesh{
  /*
   * calculate En(kx,ky)
   * fourier transform of En(kx,ky)
   * 
   */
 public:
  double px1, px2, py1, py2;
  int Lq, numOfBand;
  double eps;
  Eigen<T> *esol;  
  std::vector<std::vector<D2MeshPt > > grid;
  D2Mesh(Eigen<T> *eig, double x1, double x2, double y1, double y2, int b, int L=128,double dk=0.00001):
      esol(eig), px1(x1),px2(x2),py1(y1),py2(y2),numOfBand(b),Lq(L),eps(dk){
    // currently, Lq must be equal to Lq;
    grid.resize(Lq);
    pt.en.resize(numOfBand);
    pt.d1.resize(numOfBand);
    pt.d2.resize(numOfBand);
    ki.resize(2);
    veci.resize(esol->hamiltonian.size());
    for(int i=0; i<veci.size(); i++)
      veci[i].resize(numOfBand);
  }
  D2Mesh(std::string &infile){
    std::ifstream in(infile.c_str());
    if(!in){
    std::cerr << "[D2Mesh] error: can't open config file" << std::endl;
    std::exit(1.0);
    }
    datFromFile(in);
  }
  
  double kx_i(int index);
  double ky_i(int index);
  void cal_grid();
  void d2FourTransform(int i=0, std::ostream &out=std::cout);  // only fourier a single band
  void datToFile(std::ostream &out);
  void datFromFile(std::istream &in);
  
  DoubleVector ki;
  std::vector<std::vector<T> > veci;
  D2MeshPt pt;
};

bool justTrue(DoubleVector &p2){
  return true;
}

template<class T>
class D2Dos{
  /*
   * DOS 
   *
   */   
 public:
  D2Mesh<T>* d2m;
  DoubleMatrix xx;
  DoubleMatrix yy;
  double emin, emax, zoneArea;
  int ne, nq, numOfBand;
  double de, dq;
  double qx1, qx2, qy1, qy2;
  double eps; 
  D2Dos(D2Mesh<T>* dm, double e1=0.0, double e2=8.0, double n=1000, bool truefunc(DoubleVector &)=&justTrue, double tol=0.0000001):
      d2m(dm),emin(e1),emax(e2),ne(n),eps(tol){
    nq=d2m->Lq;
    qx1=d2m->px1; qx2=d2m->px2; qy1=d2m->py1; qy2=d2m->py2;
    dq=(qx2-qx1)/nq;
    de=(emax-emin)/ne;
    numOfBand=d2m->numOfBand;
    xx.resize(ne);
    yy.resize(ne);
    for(int i=0; i<ne; i++){
      xx[i].resize(1+numOfBand);
      yy[i].resize(1+numOfBand);
    }
    inZoneOrNot = truefunc;
    get_ZoneArea();
  }
  void dosFunc();
  void get_ZoneArea();
  bool (*inZoneOrNot)(DoubleVector &kd2);
  void hilbertTransform();
  void xyToFile(std::ostream &out);
  // double (*MatrixElement)(DoubleVector &k1, DoubleVector &k2);
  // only for correlation function
};



#endif /*_2DLATTICE_H_*/
