#include "nb.h"
#include "nb.cpp"

int main(){
  std::string struct_name="./gnuplot/d1test";
  std::ofstream zone_en(struct_name+"_zone_en.dat");
  int dim=1;
  DoubleMatrix basis(2);
  DoubleVector coupling(3);
  std::string soc="yes";
  int numOfMagnetic=2;
  double qx, qy, qz, omega;
  qx=0.5; qy=0.0; qz=0; omega=0.00;
  for(int i=0; i<basis.size(); i++)
    basis[i].resize(dim);
  basis[0][0]=1.0;
  basis[1][0]=2.0;
  
  coupling[0]=3.0; coupling[1]=1.0; coupling[2]=-2.0;

  LatticeStruct ls(basis, coupling, dim, struct_name, soc, numOfMagnetic, omega, qx, qy, qz);

  print_matrix(ls.Sx);

  double k_max=20.0;
  Adj_k ak(ls, k_max);
  ak.generate_adj_k();
  ak.print_k_set();

  Adj_k_network akn(ls);
  akn.get_adj_network(ak);
  akn.print_adj_network();

  DoubleMatrix zone(2);  
  for(int i=0; i<zone.size(); i++)
    zone[i].resize(dim);
  zone[0][0]=-0.5;
  zone[1][0]=0.5;
  
  Eigen<double> eig(&akn,128,8);  
  eig.cal_d1zone(zone, zone_en);

}
