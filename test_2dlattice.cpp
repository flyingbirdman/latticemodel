#include "2dlattice.h"
#include "2dlattice.cpp"

int main(){
  std::string struct_name="./gnuplot/2dlattice";
  std::ofstream d2grid(struct_name+"_grid.dat");
  std::ofstream d2fft(struct_name+"_fft.dat");
  std::ofstream d2dos_en(struct_name+"_dos.dat");
  int dim=2;
  DoubleMatrix basis(2);
  DoubleVector coupling(3);
  std::string soc="yes";
  int numOfMagnetic=2;
  double qx, qy, qz, omega;
  qx=0.5; qy=0.5; qz=0; omega=0.0;
  for(int i=0; i<basis.size(); i++)
    basis[i].resize(dim);
  basis[0][0]=1.0; basis[0][1]=0.0;
  basis[1][0]=0.0; basis[1][1]=1.0;
  coupling[0]=2; coupling[1]=1.0; coupling[2]=1;

  LatticeStruct ls(basis, coupling ,dim, struct_name, soc, numOfMagnetic, omega, qx, qy, qz);

  double k_max=8.0;
  Adj_k ak(ls, k_max);
  ak.generate_adj_k();

  Adj_k_network akn(ls);
  akn.get_adj_network(ak);
  
  if(soc=="yes"&& qy>1.0e-10){
    Eigen<std::complex<double> > eig(&akn,32,8);
    //eig.cal_loop(loop,loop_en);
    D2Mesh<std::complex<double> > d2m(&eig,-0.5,0.5,-0.5,0.5,8,32);
    d2m.cal_grid();
    d2m.d2FourTransform(0,d2fft);
    d2m.datToFile(d2grid);
    D2Dos<std::complex<double> > d2d(&d2m,0, 8, 1000);
    d2d.dosFunc();
    d2d.hilbertTransform();
    d2d.xyToFile(d2dos_en);
  }
  else{
    Eigen<double> eig(&akn,32,8);
    //eig.cal_loop(loop, loop_en);
    D2Mesh<double> d2m(&eig,-0.5,0.5,-0.5,0.5,8,64);
    d2m.cal_grid();
    d2m.d2FourTransform(0,d2fft);
    d2m.datToFile(d2grid);
    D2Dos<double> d2d(&d2m,0, 8, 1000);
    d2d.dosFunc();
    d2d.hilbertTransform();
    d2d.xyToFile(d2dos_en);
  }
}
