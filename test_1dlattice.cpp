#include "1dlattice.h"
#include "1dlattice.cpp"

int main(){
  std::string struct_name="./gnuplot/1dlattice";
  std::ofstream zone_en("d1test_zone_en.dat");
  std::ofstream d1grid(struct_name+"_grid.dat");
  std::ofstream d1fft(struct_name+"_fft.dat");
  std::ofstream d1dos_en(struct_name+"_dos.dat");
  int dim=1;
  DoubleMatrix basis(2);
  DoubleVector coupling(3);
  std::string soc="no";
  int numOfMagnetic=2;
  double qx, qy, qz, omega;
  qx=0.; qy=0; qz=0; omega=0.0;
  for(int i=0; i<basis.size(); i++)
    basis[i].resize(dim);
  basis[0][0]=1.0; 
  basis[1][0]=2.0; 
  coupling[0]=2; coupling[1]=1.0; coupling[2]=1.0;

  LatticeStruct ls(basis, coupling ,dim, struct_name, soc, numOfMagnetic, omega, qx, qy, qz);

  double k_max=14.0;
  Adj_k ak(ls, k_max);
  ak.generate_adj_k();

  Adj_k_network akn(ls);
  akn.get_adj_network(ak);

  DoubleMatrix zone(2);  
  for(int i=0; i<zone.size(); i++)
    zone[i].resize(dim);
  zone[0][0]=-0.5;
  zone[1][0]=0.5;

  if(soc=="yes"&& qy>1.0e-10){
    Eigen<std::complex<double> > eig(&akn,32,8);
    eig.cal_d1zone(zone,zone_en);
    D1Mesh<std::complex<double> > d1m(&eig,-0.5,0.5,8,32);
    d1m.cal_grid();
    d1m.d1FourTransform(0,d1fft);
    d1m.datToFile(d1grid);
    D1Dos<std::complex<double> > d1d(&d1m,0, 8, 2000);
    d1d.dosFunc();
    d1d.hilbertTransform();
    d1d.xyToFile(d1dos_en);
  }
  else{
    Eigen<double> eig(&akn,32,8);
    eig.cal_d1zone(zone,zone_en);
    D1Mesh<double> d1m(&eig,-0.5,0.5,8,64);
    d1m.cal_grid();
    d1m.d1FourTransform(0,d1fft);
    d1m.datToFile(d1grid);    
    D1Dos<double> d1d(&d1m,0, 8, 2000);
    d1d.dosFunc();
    d1d.hilbertTransform();
    d1d.xyToFile(d1dos_en);
  }
}
