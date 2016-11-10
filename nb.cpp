#include "nb.h"
#include "matrix_operation.h"
#include "vector_operation.h"

void LatticeStruct::generate_spin_matrix(){
  double s=(numOfMagnetic-1.0)/2.0;
  Sx.resize(numOfMagnetic); Sy.resize(numOfMagnetic); Sz.resize(numOfMagnetic);  
  for(int i=0; i<numOfMagnetic; i++){
    Sx[i].resize(numOfMagnetic); Sy[i].resize(numOfMagnetic); Sz[i].resize(numOfMagnetic);
  }
  for(int i=1; i<numOfMagnetic; i++){
    Sx[i][i-1]=sqrt(2.0*i*(s+1.0)-(i+1.0)*i)/2.0;
    Sx[i-1][i]=Sx[i][i-1];
    Sy[i][i-1]=std::complex<double>(0,sqrt(2.0*i*(s+1.0)-(i+1.0)*i)/2.0);
    Sy[i-1][i]=std::conj(Sy[i][i-1]);
    Sz[i][i]=s-i;  
    }
  Sz[0][0]=s;
}
void Adj_k::generate_adj_k(){
  double minlength=1000;
  double temp;
  int max_index;
  DoubleVector pt(ls.dim);  
  for(int i=0; i<ls.basis.size(); i++){
    temp=length(ls.basis[i]);
    if(temp<minlength)
      minlength=temp;
  }
  max_index= int(k_max/minlength);
  if(max_index <1.0e-10){
    std::cout << "[Adj_k]: one of basis is too small" << std::endl;
    std::exit(1.0);
  }
  if(ls.numOfBasis==1){
    for(int i=-max_index; i<=max_index; i++){
      for(int vi=0; vi<ls.dim; vi++)
        pt[vi]=i*ls.basis[0][vi];
      if(inRange(pt))
        k_set.insert(pt);            
    }
  }else if(ls.numOfBasis==2){
    for(int i=-max_index; i<=max_index; i++)
      for(int j=-max_index; j<=max_index; j++){
        for(int vi=0; vi<ls.dim; vi++)
          pt[vi]=i*ls.basis[0][vi]+j*ls.basis[1][vi];        
        if(inRange(pt))
          k_set.insert(pt);        
      }
  }else if(ls.numOfBasis==3){
    for(int i=-max_index; i<=max_index; i++)
      for(int j=-max_index; j<=max_index; j++)
        for(int k=-max_index; k<=max_index; k++){
          for(int vi=0; vi<ls.dim; vi++)
            pt[vi]=i*ls.basis[0][vi]+j*ls.basis[1][vi]+k*ls.basis[2][vi];
          if(inRange(pt))
            k_set.insert(pt);
        }
  }else if(ls.numOfBasis==4){
    for(int i1=-max_index; i1<=max_index; i1++)
      for(int i2=-max_index; i2<=max_index; i2++)
        for(int j1=-max_index; j1<=max_index; j1++)
          for(int j2=-max_index; j2<=max_index; j2++){
            for(int vi=0; vi<ls.dim; vi++)
              pt[vi]=i1*ls.basis[0][vi]+i2*ls.basis[1][vi]+j1*ls.basis[2][vi]+j2*ls.basis[3][vi];
            if(inRange(pt))
              k_set.insert(pt);
          }
  }else if(ls.numOfBasis==5){
    for(int i1=-max_index; i1<=max_index; i1++)
      for(int i2=-max_index; i2<=max_index; i2++)
        for(int j1=-max_index; j1<=max_index; j1++)
          for(int j2=-max_index; j2<=max_index; j2++)
            for(int k=-max_index; k<=max_index; k++){
              for(int vi=0; vi<ls.dim; vi++)
                pt[vi]=i1*ls.basis[0][vi]+i2*ls.basis[1][vi]+j1*ls.basis[2][vi]+j2*ls.basis[3][vi]
                    +k*ls.basis[4][vi];
              if(inRange(pt))
                k_set.insert(pt);
          }
  }else if(ls.numOfBasis==6){
    for(int i1=-max_index; i1<=max_index; i1++)
      for(int i2=-max_index; i2<=max_index; i2++)
        for(int j1=-max_index; j1<=max_index; j1++)
          for(int j2=-max_index; j2<=max_index; j2++)
            for(int k1=-max_index; k1<=max_index; k1++)
              for(int k2=-max_index; k2<=max_index; k2++){
                for(int vi=0; vi<ls.dim; vi++)
                  pt[vi]=i1*ls.basis[0][vi]+i2*ls.basis[1][vi]+j1*ls.basis[2][vi]
                      +j2*ls.basis[3][vi]+k1*ls.basis[4][vi]+k2*ls.basis[5][vi];
              if(inRange(pt))
                k_set.insert(pt);
          }
  }else if(ls.numOfBasis==7){
    for(int i1=-max_index; i1<=max_index; i1++)
      for(int i2=-max_index; i2<=max_index; i2++)
        for(int j1=-max_index; j1<=max_index; j1++)
          for(int j2=-max_index; j2<=max_index; j2++)
            for(int k1=-max_index; k1<=max_index; k1++)
              for(int k2=-max_index; k2<=max_index; k2++)
                for(int l1=-max_index; l1<=max_index; l1++){
                  for(int vi=0; vi<ls.dim; vi++)
                    pt[vi]=i1*ls.basis[0][vi]+i2*ls.basis[1][vi]+j1*ls.basis[2][vi]
                        +j2*ls.basis[3][vi]+k1*ls.basis[4][vi]+k2*ls.basis[5][vi]
                        +l1*ls.basis[6][vi];
              if(inRange(pt))
                k_set.insert(pt);
          }
  }else if(ls.numOfBasis==8){
    for(int i1=-max_index; i1<=max_index; i1++)
      for(int i2=-max_index; i2<=max_index; i2++)
        for(int j1=-max_index; j1<=max_index; j1++)
          for(int j2=-max_index; j2<=max_index; j2++)
            for(int k1=-max_index; k1<=max_index; k1++)
              for(int k2=-max_index; k2<=max_index; k2++)
                for(int l1=-max_index; l1<=max_index; l1++)
                  for(int l2=-max_index; l2<=max_index; l2++){
                    for(int vi=0; vi<ls.dim; vi++)
                      pt[vi]=i1*ls.basis[0][vi]+i2*ls.basis[1][vi]+j1*ls.basis[2][vi]
                          +j2*ls.basis[3][vi]+k1*ls.basis[4][vi]+k2*ls.basis[5][vi]
                          +l1*ls.basis[6][vi]+l2*ls.basis[7][vi];
                    if(inRange(pt))
                      k_set.insert(pt);
                  }
  }else{
    std::cerr << "[Adj_k]: too many lattice basis, max # of basis is 8. You need to recode" << std::endl;
    std::exit(1);
  }
  if(k_set.size()<1){
    std::cerr << "[Adj_k] generate_adj_k:" << k_set.size() << " something in lattice config must be wrong! " << std::endl;
    std::exit(1.0);
  }
  if(ls.socOrNot=="yes" || ls.socOrNot=="y")
    dimOfHilbert=k_set.size()*ls.numOfMagnetic;
  else
    dimOfHilbert=k_set.size();
}

void Adj_k::print_k_set(){
  std::clog << "# " << ls.struct_name << ":  all the coupling points in k space :#" << std::endl;
  cmpSet::iterator iter=k_set.begin();
  for(;iter != k_set.end(); ++iter){
    for(int i=0; i < ls.dim; ++i)
      std::clog << (*iter)[i] << "  ";
    std::clog << std::endl;
  }
  std::clog << "# ********************************************************** #" << std::endl;
}

void Adj_k_network::get_adj_network(Adj_k &ak){
  cmpSet::iterator iter=ak.k_set.begin();
  cmpSet::iterator iter_find;
  DoubleVector adj_vec(dim); 
  std::list<std::pair<int,double> > adj_indx;
  int index;

  for(int i=0; iter !=ak.k_set.end(); ++i, ++iter){
    DoubleVector cur_vec(*iter);
    adj_indx.push_back(std::pair<int,double>(i,ls.coupling[0]));
    for(int j=0; j < ls.numOfBasis;j++){ 
      for(int k=0; k < dim; k++)
        adj_vec[k]=cur_vec[k]+ls.basis[j][k]; //case one: +basis
      iter_find=(ak.k_set).find(adj_vec);
      if(iter_find != ak.k_set.end()){
        index=std::distance(ak.k_set.begin(), iter_find);        
        adj_indx.push_back(std::pair<int,double> (index, ls.coupling[j+1]));
      }
      else{
        if(*iter==adj_vec){
          index=std::distance(ak.k_set.begin(), iter_find);
          adj_indx.push_back(std::pair<int,double> (index, ls.coupling[j+1]));
        }
      }
    }
    
    for(int j=0; j < ls.numOfBasis;j++){ 
      for(int k=0; k < dim; k++)
        adj_vec[k]=cur_vec[k]-ls.basis[j][k]; //case two: -basis
      iter_find=(ak.k_set).find(adj_vec);
      if(iter_find != (ak.k_set).end()){
        index=std::distance(ak.k_set.begin(), iter_find);        
        adj_indx.push_back(std::pair<int,double> (index, ls.coupling[j+1]));
      }
      else{
        if(*iter==adj_vec){
          index=std::distance(ak.k_set.begin(), iter_find);
          adj_indx.push_back(std::pair<int,double> (index, ls.coupling[j+1]));
        }
      }
    }
    adj_network.push_back(Adj_list<double> (cur_vec, adj_indx));
    adj_indx.clear();
  }
  dimOfHilbert=ak.dimOfHilbert;
}

void Adj_k_network::print_adj_network(){
  std::list<std::pair<int,double> >::iterator iter;
  std::clog << "# adj_network of the structure :" << ls.struct_name << std::endl; 
  for(int i=0; i<adj_network.size(); i++){
       for(iter=adj_network[i].indx.begin(); iter != adj_network[i].indx.end(); ++iter)
            std::clog << "[" << (*iter).first << ", " << (*iter).second << "] "; 
       std::clog << std::endl;
  }
  std::clog << "# ********************************************************** #" << std::endl;
}

template<class T>
void Eigen<T>::cal_ei(DoubleVector &ki, DoubleVector &val, std::vector<std::vector<T> > &vec, char jobz){
  get_diag_hamiltonian(ki, hamiltonian, *akn);
  eigen_symmetric_solver(hamiltonian, val, vec,jobz);
}

template<class T>
void Eigen<T>::cal_loop(DoubleMatrix &loop, std::ofstream &out){
  DoubleVector k12(dim);
  DoubleVector kij(dim);
  double length;
  double totallength = 0;
  double step = 0;

  for(int i=0; i < loop.size()-1; i++){
    length = 0;
    for(int j=0; j < dim; j++){
      k12[j]=loop[i+1][j]-loop[i][j];
      length+= k12[j]*k12[j];
    }
    totallength += sqrt(length);
  }

  for(int i=0; i < loop.size()-1; i++){
    length=0;
    for(int j=0; j < dim; j++){
      k12[j]=loop[i+1][j]-loop[i][j];
      length+= k12[j]*k12[j];
    }
    length=sqrt(length);
    step+=0.5*length/totallength/fftsize;
    for(int j=0; j < fftsize; j++){
      for(int k=0; k < dim; k++)
        kij[k] = loop[i][k]+k12[k] * (j+0.5) / fftsize;       
      out << step << "  ";
      get_diag_hamiltonian(kij, hamiltonian, *akn);
      eigen_symmetric_solver(hamiltonian,eigval,eigvec);
      obtain_eigval_dat(kij,out);
      step += length/totallength/fftsize;
    }
    step-=0.5*length/totallength/fftsize;
  }
}

template<class T>
void Eigen<T>::cal_d1zone(DoubleMatrix &zone, std::ofstream &out){
  if(zone.size()!=2 && dim!=1){
    std::clog << "cal_d1zone: dim must be 1 and zone.size()=2" << std::endl;    
    std::exit(1);
  }
  double k12;
  DoubleVector kij(dim);  
  
  k12=zone[1][0]-zone[0][0];
  for(int j=0; j < fftsize; j++){
    kij[0] = zone[0][0]+ k12 * (j+0.5) / fftsize;       
    get_diag_hamiltonian(kij, hamiltonian, *akn);
    eigen_symmetric_solver(hamiltonian,eigval,eigvec);
    obtain_eigval_dat(kij,out);
  }
}

template<class T>
void Eigen<T>::cal_d2zone(DoubleMatrix &zone, std::ofstream &out){
  if(dim!=2 && zone.size()!=3){
    std::clog << "cal_d2zone: dim must be 2 and zone.size()=3" << std::endl;
    std::exit(1);
  }
  DoubleVector kx(dim);
  DoubleVector ky(dim);
  DoubleVector kij(dim);

  for(int i=0; i < dim; i++){
    kx[i]=zone[1][i]-zone[0][i];
    ky[i]=zone[2][i]-zone[0][i];
  }
  for(int i=0; i < fftsize; i++){
    for(int j=0; j < fftsize; j++){      
      for(int k=0; k < dim; k++)
        kij[k] = zone[0][k]+ kx[k]*(i+0.5)/fftsize + ky[k]*(j+0.5)/fftsize;
      get_diag_hamiltonian(kij, hamiltonian, *akn);
      eigen_symmetric_solver(hamiltonian,eigval,eigvec);
      obtain_eigval_dat(kij,out);
    }
    out << std::endl;
  }
}

template<class T>
void Eigen<T>::cal_d3zone(DoubleMatrix &zone, std::ofstream &out){
  if(dim!=3 && zone.size()!=4){
    std::clog << "cal_d3zone: dim must be 3 and zone.size()=4" << std::endl;
    std::exit(1);
  }
  DoubleVector kx(dim);
  DoubleVector ky(dim);
  DoubleVector kz(dim);
  DoubleVector kij(dim);

  for(int i=0; i < dim; i++){
    kx[i]=zone[1][i]-zone[0][i];
    ky[i]=zone[2][i]-zone[0][i];
    kz[i]=zone[3][i]-zone[0][i];
  }
  for(int i=0; i < fftsize; i++){
    for(int j=0; j < fftsize; j++){
      for(int k=0; k < fftsize; k++){      
        for(int vi=0; vi < dim; vi++)
          kij[vi] = zone[0][vi]+ kx[vi]*(i+0.5)/fftsize + ky[vi]*(j+0.5)/fftsize+kz[vi]*(k+0.5)/fftsize;
        get_diag_hamiltonian(kij, hamiltonian, *akn);
        eigen_symmetric_solver(hamiltonian,eigval,eigvec);
        obtain_eigval_dat(kij,out);
      }
      out << std::endl;
    }
    out << std::endl << std::endl;
  }
}

template<class T>
void Eigen<T>::obtain_eigval_dat(const DoubleVector &k, std::ofstream &out){
  for(int i=0; i < dim; i++)
    out << k[i] << "  ";
  for(int i=0; i < eigval.size(); i++)
    out << eigval[i] << "  ";
  out << std::endl;
}

template <class T>
void SubEigen<T>::cal_lk_subspace(DoubleVector &ki){  
  get_diag_hamiltonian(ki, h1, *akn1);
  get_diag_hamiltonian(ki, h2, *akn2);
  eigen_symmetric_solver(h1, val1, vec1,'V',bl);
  eigen_symmetric_solver(h2, val2, vec2,'V',bl);
}

template <class T>
void SubEigen<T>::cal_d1_subspace(DoubleMatrix &zone, std::ofstream &out){
  if(zone.size()!=2 && dim!=1){
    std::clog << "cal_d1_subspace: dim must be 1 and zone.size()=2" << std::endl;    
    std::exit(1);
  }
  double k12;
  DoubleVector kij(dim);  
  
  k12=zone[1][0]-zone[0][0];
  for(int j=0; j < fftsize; j++){
    kij[0] = zone[0][0]+ k12 * (j+0.5) / fftsize; 
    cal_lk_subspace(kij);
    colms_inner_product(lvec,vec1,vec2);
    colms_inner_product_with_matrix(lh,vec1,h2,vec1);
  }
}

template <class T>
void SubEigen<T>::cal_d2_subspace(DoubleMatrix &zone, std::ofstream &out){
  if(dim!=2 && zone.size()!=3){
    std::clog << "cal_d2_subspace: dim must be 2 and zone.size()=3" << std::endl;
    std::exit(1);
  }
  DoubleVector kx(dim);
  DoubleVector ky(dim);
  DoubleVector kij(dim);

  for(int i=0; i < dim; i++){
    kx[i]=zone[1][i]-zone[0][i];
    ky[i]=zone[2][i]-zone[0][i];
  }
  for(int i=0; i < fftsize; i++){
    for(int j=0; j < fftsize; j++){      
      for(int k=0; k < dim; k++)
        kij[k] = zone[0][k]+ kx[k]*(i+0.5)/fftsize + ky[k]*(j+0.5)/fftsize;
      cal_lk_subspace(kij);
      colms_inner_product(lvec,vec1,vec2);
      colms_inner_product_with_matrix(lh,vec1,h2,vec1);
    }
  }
}

template <class T>
void SubEigen<T>::cal_d3_subspace(DoubleMatrix &zone, std::ofstream &out){
  if(dim!=3 && zone.size()!=4){
    std::clog << "cal_d3_subspace: dim must be 3 and zone.size()=4" << std::endl;
    std::exit(1);
  }
  DoubleVector kx(dim);
  DoubleVector ky(dim);
  DoubleVector kz(dim);
  DoubleVector kij(dim);

  for(int i=0; i < dim; i++){
    kx[i]=zone[1][i]-zone[0][i];
    ky[i]=zone[2][i]-zone[0][i];
    kz[i]=zone[3][i]-zone[0][i];
  }
  for(int i=0; i < fftsize; i++){
    for(int j=0; j < fftsize; j++){
      for(int k=0; k < fftsize; k++){      
        for(int vi=0; vi < dim; vi++)
          kij[vi] = zone[0][vi]+ kx[vi]*(i+0.5)/fftsize + ky[vi]*(j+0.5)/fftsize+kz[vi]*(k+0.5)/fftsize;
        cal_lk_subspace(kij);
        colms_inner_product(lvec,vec1,vec2);
        colms_inner_product_with_matrix(lh,vec1,h2,vec1);
      }
    }
  }
}

void get_off_diag_hamiltonian(DoubleMatrix &h, const  Adj_k_network &akn){
  std::list<std::pair<int,double> >::const_iterator iter;
  if(akn.ls.socOrNot=="yes"||akn.ls.socOrNot=="y"){
    int m=akn.ls.numOfMagnetic;
    for(int i=0; i < akn.adj_network.size(); i++)
      for(iter =akn.adj_network[i].indx.begin(); iter != akn.adj_network[i].indx.end(); ++iter)
        for(int j=0; j<m; j++){
          h[i*m+j][(iter->first)*m+j] = iter->second * 0.5;
        }
  }else{
    for(int i=0; i < akn.adj_network.size(); i++)
      for(iter =akn.adj_network[i].indx.begin(); iter != akn.adj_network[i].indx.end(); ++iter)
        h[i][iter->first] = iter->second * 0.5;   
  }
}
void get_off_diag_hamiltonian(ComplexMatrix &h, const  Adj_k_network &akn){
  std::list<std::pair<int,double> >::const_iterator iter;  
  if(akn.ls.socOrNot=="yes"||akn.ls.socOrNot=="y"){
    int m=akn.ls.numOfMagnetic;
    for(int i=0; i < akn.adj_network.size(); i++)
      for(iter =akn.adj_network[i].indx.begin(); iter != akn.adj_network[i].indx.end(); ++iter)
        for(int j=0; j<m; j++){
          h[i*m+j][(iter->first)*m+j] = iter->second * 0.5;
        }
  }else{
    for(int i=0; i < akn.adj_network.size(); i++)
      for(iter =akn.adj_network[i].indx.begin(); iter != akn.adj_network[i].indx.end(); ++iter)
        h[i][iter->first] = iter->second * 0.5;   
  }
}
void get_diag_hamiltonian(const DoubleVector &p, DoubleMatrix &h, const Adj_k_network &akn){
  int n=akn.adj_network.size();
  if(akn.ls.socOrNot=="yes"||akn.ls.socOrNot=="y"){
    double Qx, Qy, Qz, Omega, e0;
    DoubleMatrix Sx(akn.ls.Sx);
    DoubleMatrix Sz(akn.ls.Sz);
    DoubleMatrix Sx2(Sx);
    int m=akn.ls.numOfMagnetic;    
    for(int i=0; i<m; i++)
      for(int j=0; j<m; j++){
        Sx2[i][j]=0;
        for(int k=0; k<m; k++)
          Sx2[i][j]+=Sx[i][k]*Sx[k][j];
      }
    Qx=akn.ls.Qx; Qy=akn.ls.Qy; Qz=akn.ls.Qz; Omega=akn.ls.Omega;
    if(fabs(Qy)>1.0e-10){
      std::cerr << "[double] get_diag_hamiltonian: using wrong data type" << std::endl;
      std::exit(1.0);
    }
    e0=akn.ls.coupling[0];
    if(fabs(Qx)>1.0e-10)
      for(int i=0; i < n; i++)
        for(int j1=0; j1 < m; j1++)
          for(int j2=0; j2 < m; j2++)
            h[i*m+j1][i*m+j2] = -2.0*Qx*Sx[j1][j2]*(p[0]+akn.adj_network[i].vec[0])+Qx*Qx*Sx2[j1][j2];
    else
      for(int i=0; i < n; i++)
        for(int j1=0; j1 < m; j1++)
          for(int j2=0; j2 < m; j2++)
            h[i*m+j1][i*m+j2] = 0;
    if(akn.dim==3 && fabs(Qz)>1.0e-10)
      for(int i=0; i < n; i++)
        for(int j=0; j < m; j++)
          h[i*m+j][i*m+j] += -2.0*Qz*Sz[j][j]*(p[2]+akn.adj_network[i].vec[2])+Qz*Qz*Sz[j][j]*Sz[j][j]
              ;
    for(int i=0; i < n; i++)
      for(int j=0; j<m; j++){
        h[i*m+j][i*m+j] += e0+Omega*Sz[j][j];
        for(int k=0; k < akn.dim; k++)
          h[i*m+j][i*m+j] += (p[k]+akn.adj_network[i].vec[k])*(p[k]+akn.adj_network[i].vec[k]);
      }
  }else{
    double e0=akn.ls.coupling[0];
    for(int i=0; i < n; i++){
      h[i][i] = e0;
      for(int j=0; j < akn.dim; j++)
        h[i][i] += (p[j]+akn.adj_network[i].vec[j])*(p[j]+akn.adj_network[i].vec[j]);
    }
  }
}

void get_diag_hamiltonian(const DoubleVector &p, ComplexMatrix &h, const Adj_k_network &akn){
  int n=akn.adj_network.size();
  if(akn.ls.socOrNot=="yes"||akn.ls.socOrNot=="y"){
    double Qx, Qy, Qz, Omega, e0;
    DoubleMatrix Sx(akn.ls.Sx);
    ComplexMatrix Sy(akn.ls.Sy);
    DoubleMatrix Sz(akn.ls.Sz);
    DoubleMatrix Sx2(Sx);
    ComplexMatrix Sy2(Sy);
    int m=akn.ls.numOfMagnetic;    
    for(int i=0; i<m; i++)
      for(int j=0; j<m; j++){
        Sx2[i][j]=0;
        Sy2[i][j]=0;
        for(int k=0; k<m; k++){
          Sx2[i][j]+=Sx[k][i]*Sx[k][j];
          Sy2[i][j]+=Sy[i][k]*Sy[k][j];
        }
      }
    Qx=akn.ls.Qx; Qy=akn.ls.Qy; Qz=akn.ls.Qz; Omega=akn.ls.Omega; e0=akn.ls.coupling[0];
    if(fabs(Qx)>1.0e-10)
      for(int i=0; i < n; i++)
        for(int j1=0; j1 < m; j1++)
          for(int j2=0; j2 < m; j2++)
            h[i*m+j1][i*m+j2] = -2.0*Qx*Sx[j1][j2]*(p[0]+akn.adj_network[i].vec[0])+Qx*Qx*Sx2[j1][j2];
    else
      for(int i=0; i < n; i++)
        for(int j1=0; j1 < m; j1++)
          for(int j2=0; j2 < m; j2++)
            h[i*m+j1][i*m+j2] = 0;
    if(akn.dim>=2 && fabs(Qy)>1.0e-10)
      for(int i=0; i < n; i++)
        for(int j1=0; j1 < m; j1++)
          for(int j2=0; j2 < m; j2++)
            h[i*m+j1][i*m+j2] += -2.0*Qy*Sy[j1][j2]*(p[1]+akn.adj_network[i].vec[1])+Qy*Qy*Sy2[j1][j2];
    if(akn.dim==3 && fabs(Qz)>1.0e-10)
      for(int i=0; i < n; i++)
        for(int j=0; j < m; j++)
          h[i*m+j][i*m+j] += -2.0*Qz*Sz[j][j]*(p[2]+akn.adj_network[i].vec[2])+Qz*Qz*Sz[j][j]*Sz[j][j];
    for(int i=0; i < n; i++)
      for(int j=0; j<m; j++){
        h[i*m+j][i*m+j] += e0+Omega*Sz[j][j];;
        for(int k=0; k < akn.dim; k++)
          h[i*m+j][i*m+j] += (p[k]+akn.adj_network[i].vec[k])*(p[k]+akn.adj_network[i].vec[k]);
      }
  }else{
    double e0=akn.ls.coupling[0];
    for(int i=0; i < n; i++){
      h[i][i] = e0;
      for(int j=0; j < akn.dim; j++)
        h[i][i] += (p[j]+akn.adj_network[i].vec[j])*(p[j]+akn.adj_network[i].vec[j]);
    }
  }
}
