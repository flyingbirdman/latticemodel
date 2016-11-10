#include "nb.h"
#include "nb.cpp"

int main(){
  std::string struct_name="./gnuplot/d2test";
  std::ofstream loop_en(struct_name+"_loop_en.dat");
  std::ofstream zone_en(struct_name+"_zone_en.dat");
  int dim=2;
  DoubleMatrix basis(2);
  DoubleVector coupling(3);
  std::string soc="yes";
  int numOfMagnetic=2;
  double qx, qy, qz, omega;
  qx=0.5; qy=0.5; qz=0; omega=0.0;
  for(int i=0; i<basis.size(); i++)
    basis[i].resize(dim);
  basis[0][0]=1.0; basis[0][1]=0;
  basis[1][0]=0.0; basis[1][1]=1.0;
  
  coupling[0]=2; coupling[1]=1.0; coupling[2]=1.;

  LatticeStruct ls(basis, coupling, dim, struct_name, soc, numOfMagnetic, omega, qx, qy, qz);

  double k_max=8.0;
  Adj_k ak(ls, k_max);
  ak.generate_adj_k();
  ak.print_k_set();

  Adj_k_network akn(ls);
  akn.get_adj_network(ak);
  akn.print_adj_network();

  DoubleMatrix loop(4);
  DoubleMatrix zone(3);
  for(int i=0; i<loop.size(); i++)
    loop[i].resize(dim);
  for(int i=0; i<zone.size(); i++)
    zone[i].resize(dim);
  loop[0][0]=0.0; loop[0][1]=0.0;
  loop[1][0]=0.5; loop[1][1]=0.0;
  loop[2][0]=0.5; loop[2][1]=0.5;
  loop[3][0]=0.0; loop[3][1]=0.0;
  zone[0][0]=-0.5; zone[0][1]=-0.5;
  zone[1][0]=0.5; zone[1][1]=-0.5;
  zone[2][0]=-0.5; zone[2][1]=0.5;
  
  if(soc=="yes"&& qy>1.0e-10){
    Eigen<std::complex<double> > eig(&akn,64,8);
    std::cout << "complex" << std::endl;
    eig.cal_loop(loop,loop_en);
    eig.cal_d2zone(zone, zone_en);
  }
  else{
    Eigen<double> eig(&akn,64,8);
    eig.cal_loop(loop, loop_en);
    eig.cal_d2zone(zone, zone_en);
  }
}
